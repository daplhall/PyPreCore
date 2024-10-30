#%matplotlib widget
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
import scipy.optimize as sciopt
import scipy.interpolate as scilate
import seaborn as sns
sns.set(context = 'notebook',font = 'serif',font_scale = 1.3)
sns.set_palette('Set2', 8)


import scaling as sc
import os
def magnitude(r, keepdims = False):
    """ 
    scalar magnitude of a vector
    from num astro
    """
    return np.sqrt(np.sum(r**2,1,keepdims=keepdims))

def cube_function(func):
    def function_cubed(r):
        return func(r)**3
    return function_cubed

class Particles(object):
    __slots__ = 'n','t', 'dt' , 'env','units','r','s', 'Stmin' , 'Stmax'\
                ,'rho_d', 'v','M' ,'D', 'ngrains','sinterface', 'runmask'
    def __init__(self, n, env, R, Stmin, Stmax,rho_d = 1.6, C = 0.004):
        self.n      = n
        self.t      = np.zeros((n     ,1),dtype = np.float64)
        self.dt     = np.zeros((self.n,1),dtype = np.float64)
        self.env    = env ## enviorment ie sphere.
        self.units  = env.units
        self.r      = np.zeros((n,2),dtype = np.float64)
        self.v      = self.r.copy()
        self.r[:,0] = R
        self.rho_d  = rho_d # g/cm ^3
        st     = np.geomspace((Stmin,), (Stmax,), n)
        rmag   = magnitude(self.r, True)
        tff    = np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rmag)))                           # TODO move to enviorment
        vth    = np.sqrt(8/np.pi * self.units.k_B*self.env.T(rmag)/(self.env.mu*self.units.m_p)) # TODO move to enviorment instead, and then just call it.
        self.s = tff * self.env.rho(rmag)* vth / self.rho_d * st
        #to ensure that we only evolve the ones that aren't done
        self.runmask = np.ones(self.n, dtype = bool)

    @property
    def search(self):
        rm = magnitude(self.r)
        return np.int64(self.env.interpolatesearch[self.env.idstep](rm))


    @property
    def St(self):
        rmag = magnitude(self.r, True)
        ts   = 1/self.env.drag(self.r,self.s,self.rho_d)
        tff  = np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rmag))) 
        return ts/tff


    def courant(self, Cdr = 0.2, Cff = 0.1, Csnap = 0.1, Ckep = 0.1, same_dt = False, verbose = False):
        
        ### Small hack to make dr work for cells in ncells -1 as it width should be inf ####
            ### This is the lazy way, the better way would be append to env._dr, but the i have to change  ###
            ### more in the genreal code base                                                              ##
        dr          = np.ones_like(self.t)*np.inf
        search      = self.search
        smask       = (search != (self.env.ncells - 1))
        dr[smask,0] = self.env.dr[search[smask]]
        dr = dr[self.runmask]
        ############################
        #dr   = self.env.dr[self.search,None][self.runmask] # prolematic if we are at inf
        rm   = magnitude(self.r[self.runmask], True)
        a    = self.gravity() #- self.env.drag(self.r,self.s,self.rho_d)*(self.v - self.env.vg(self.r))# NOTE   g - c_d(vd - vg)
        dt_dr                 = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        if np.any(self.env.rho(rm) == 0):
            debug = True
        dt_ff                 = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap               = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep                = Ckep  * np.sqrt(rm*rm*rm/(self.units.G*self.env.M(rm)))       # kepler based
        dt1                   = np.minimum(dt_dr, dt_ff) 
        dt2                   = np.minimum(dt_snap,dt_kep)
        self.dt[self.runmask] =  np.minimum(dt1, dt2)
        if verbose:
            imin = np.argmin( self.dt )
            print( 'idx: %i,  radius: %.3e, cell: %i, dt_dr: %.3e, dt_ff: %.3e, dt_snap: %.3e, dt_kep: %.3e, size: %.3e, a_grav: %.3e' \
                   %( imin, magnitude(self.r)[imin]/self.units.AU, self.search[imin], dt_dr[imin]/self.units.yr,
                   dt_ff[imin]/self.units.yr, dt_snap/self.units.yr, dt_kep[imin]/self.units.yr , self.s[imin], a[imin][0]) )
        if same_dt:
            self.dt[:] = self.dt.min()
        

    def gravity(self):
        rm = magnitude(self.r, True)
        return  - self.r * self.units.G*self.env.M(rm)/(rm*rm*rm) 
        #return  - self.r * self.units.G*self.units.m_Sun/(rm*rm*rm) # TODO add mass change
        #return -self.r*self.G/(rm*rm*rm) # self .G is temp

    def kick(self,step=0.5):
        """
        Must be paired with self.v

        returns a step for v based on gravity        
        
        """
        force  = self.gravity()#not really a force, more so an acceleration
        self.v = self.v + step*self.dt*force
        #return             step*self.dt*force
    
    def drift(self,step=1.0):
        """ 
        Must be paired with self.v

        returns a drift step for r based on velocity v.
        """
        self.r = self.r + step*self.dt*self.v
        #return             step*self.dt*self.v

    def KDK(self,same_dt=True):
        self.kick(0.5)
        self.drift(1.0)
        self.kick(0.5)

    def KDK_drag(self,same_dt=True):
        """
        TODO Important!! Rewrite KDK as it cant make 
        self.v = self.v  - Cdt*(self.v - self.env.vg(self.r))/(1 + Cdt)
        fit
        self.v = self.v  + (dtg - dt Cd(vd-vg) )/(1+cdt)
        """
        ###
        self.kick(0.5)
        # NOTE + != * :(
        r      = self.r  # + 0.5 * self.dt* self.v # half here as it needs to sync with kick
        Cdt    = 0.5*self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
        self.v = self.v  - Cdt*(self.v - self.env.vg( r ))/(1 + Cdt) # vg = 0
        ###
        self.drift(1.0)
        ###
        self.kick(0.5)
        Cdt    = 0.5*self.dt * self.env.drag(self.r,self.s)
        self.v = self.v  - Cdt*(self.v - self.env.vg(self.r))/(1 + Cdt) # vg = 0


    def DKD_drag(self,same_dt=True):
        """
        TODO Important!! Rewrite KDK as it cant make 
        self.v = self.v  - Cdt*(self.v - self.env.vg(self.r))/(1 + Cdt)
        fit
        self.v = self.v  + (dtg - dt Cd(vd-vg) )/(1+cdt)
        """
        self.drift(0.5)
        # NOTE + != * :(
        r      = self.r  # + 0.5 * self.dt* self.v # half here as it needs to sync with kick
        Cdt    = self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
        self.v = self.v  - Cdt*(self.v - self.env.vg( r ))/(1 + Cdt) # vg = 0
        self.kick(1.0)
        self.drift(0.5)



    def KDK_drag_(self):
        dt = 0.5*self.dt
        # kick
        g  = self.gravity()
        r  = self.r  + dt* self.v
        cdt = dt * self.env.KDK_drag(r, self.s, self.rho_d)
        self.v = self.v + ( dt*g - cdt * (self.v - self.env.rho( r )))/(1 + cdt)
        # drift
        self.r = self.r + 1.0*self.dt*self.v
        #kick
        g  = self.gravity()
        Cd = self.env.KDK_drag(self.r, self.s, self.rho_d)
        self.v = self.v + ( dt*g - cdt * (self.v - self.env.rho( self.r )))/(1 + cdt)

class DataSphere(object):
    """ # rc, vg, rhog , M, T 
        #  0   1     2   3  4  """
        # NOTE only works on .npy files as of now (hardcoded)
        # TODO add time interpoaltion
    __slots__ = '_rc', '_vg' ,'_rho', '_M' , '_T','_re', 'nsteps', 'ncells', 'idstep', 'dt', 't', 'interpolatevg',\
         'interpolaterho', 'interpolateM', 'interpolateT', 'mu','units' , 'interpolatesearch', '_dr','_cs','interpolateCS'
    def __init__(self, loadpath, units, idstep = 0, mu = 2.42):
        time         = np.load(loadpath + 'time.npy')
        data         = np.load(loadpath + 'data.npy')
        shape        = data.shape
        self.nsteps  = shape[-1]
        self.ncells  = shape[ 0]
        self._rc     = data[:,0,:]
        self._vg     = data[:,1,:]
        self._rho    = data[:,2,:]
        self._M      = data[:,3,:]
        self._T      = data[:,4,:]
        self._re     = data[:,5,:]
        self._cs     = data[:,6,:]
        self._dr     = data[1:,5,:] - data[:-1,5,:]
        ######
        #infi = np.ones(814)*np.inf ## NOTE TEST FOR IF PARTICLE IS IN THE LAST CELL
        #self._dr = np.vstack((self._dr,infi))
        #######
        self.idstep  = idstep # what sna
        self.dt      = time[1]
        self.t       = time[0]
        self.mu      = mu
        self.units   = units
        # interpolation arrays
        self.update()

    def update(self):
        n = range(self.ncells)                                                                         # (self._vg [0,idx],              0)
        ### vg exterpolate, what happens if r < rc[-1]?
        self.interpolatevg     = [scilate.interp1d(self._rc[:,idx], self._vg[:,idx]        ,fill_value = 'extrapolate'                                       , bounds_error = False)                  for idx in range(self.nsteps)]
        self.interpolaterho    = [scilate.interp1d(self._rc[:,idx], self._rho[:,idx]       ,fill_value = (self._rho[0,idx],              0)                  , bounds_error = False)                  for idx in range(self.nsteps)]
        interpolateM           = [scilate.interp1d(self._re[:,idx], np.cbrt(self._M[:,idx]),fill_value = (np.cbrt(self._M  [0,idx]),np.cbrt(self._M[-1,idx])), bounds_error = False, kind = 'linear') for idx in range(self.nsteps)] # re instead of rc i misunderstoond
        self.interpolateM      = [cube_function(func) for func in interpolateM]            
        self.interpolateT      = [scilate.interp1d(self._rc[:,idx], self._T[:,idx]         ,fill_value = (self._T  [0,idx],              0)                  , bounds_error = False)                  for idx in range(self.nsteps)]
        self.interpolatesearch = [scilate.interp1d(self._re[:,idx], n                      ,fill_value = (0               ,self.ncells - 1)                  , bounds_error = False)                  for idx in range(self.nsteps)]
        self.interpolateCS     = [scilate.interp1d(self._rc[:,idx], self._cs[:,idx]        ,fill_value = (self._cs [0,idx],              0)                  , bounds_error = False)                  for idx in range(self.nsteps)]


    @property
    def next(self):
        self.idstep = self.idstep + 1 if self.nsteps - 1 > self.idstep else self.nsteps -1
    @property
    def prev(self):
        self.idstep = self.idstep - 1 if self.idstep > 0 else 0
    @property
    def reset(self):
        self.idstep = 0
    @property
    def rc(self):
        return self._rc[:,self.idstep] 
    @property
    def re(self):
        return self._re[:,self.idstep] 
    @property
    def dr(self):
        return self._dr[:,self.idstep]

    @property
    def dtsnap(self):
        return self.dt[self.idstep]

    def vgmag(self,rmag):
         return np.abs(self.interpolatevg[self.idstep](rmag))

    def vg(self, r):
        """
        TODO NOTE this needs to be a 2d vector, in cartesian corodiantes, this is the radial
        """
        rm = magnitude(r,True)
        vm = self.interpolatevg[self.idstep](rm) # NOTE this is signed
        #rm[rm == 0] = 1 # TODO NOTE if rm is 0 then r is zero everywhere and thus v is 0
        return vm * r/rm

    def rho(self, rmag):
        return self.interpolaterho[self.idstep](rmag)
    
    def M(self, rmag):
        return self.interpolateM[self.idstep](rmag)
    
    def T(self, rmag):
        return self.interpolateT[self.idstep](rmag)

    def Cs(self,rmag):
        return self.interpolateCS[self.idstep](rmag)

    # Copied from sphere class
    def drag(self, r, s,rho_d = 1.6):
        """
        This is actually the inverse stopping time.
        so drag = 1/t_s
        rho_d = 1.6 g/cm^3
        """
        rm    = magnitude(r,True)
        rho_g = self.rho(rm)
        vth   = np.sqrt(8/np.pi)*self.Cs(rm)
        return rho_g*vth/(rho_d*s)

def GrainDistIntegrated(a1, a2):
    """
    This here is what the distributen  dn/da = N * n^-3.5 => dn/dlna = N * n^-2.5

    a1 is lower,min
    a2 is upper,max
    """
    return np.sqrt(a2) - np.sqrt(a1)

def Particles_grain_radial(sphere:DataSphere, func_dist, grain_interface,D = 0.01,debug = False, **kwards):
    """
    All starts along the positive x directions
    
    dust velocity is currently always 0 so same as gas to start.

    NOTE the last cell in M is M[-1] = M[ncells - 1] is not added or used, as it doens't have a n+1 over it ie a cell next to it.


    args::
        sphere          : DataSphere
        -- The sphere which contrains interpolated data from a collapse
        func_dist       : callable
        -- Solution to the integral g(a1,a2) =  int_a1^a2 dn/dlna a^2 da
           used in fraction fj = g(a_l,a_u)/g(a_min,a_amx)
        grain_interface : 1D ndarray 
        -- array of the interfaces between grain bins 
        D               : float
        -- The dust-to-gas ratio, default = 1/100
           This defualt will give truncation error as 1/100 in binary is an infinite series.
    kwards::
        rho_d : float
        -- material density of the dust, default = 1.6 g/cm^3
        C     : float
        -- Courant number, default = 0.004
    
    return: Particles Class
    -- The paricles class with grains created from a distribution
       The data is structure so [{a[:],radial bin 0},{a[:],radial bin 1},{a[:],radial bin 2},...] ie grains in radial bin 0 first then all grains in radial bin 2 etc.

    TODO BUG there seems to be a bug in placement of particles at the lower end, there are decimal diffrence 
    """
    Ngrain    = len(grain_interface) - 1 # number of grain bins not interfaces
    Nbins     = sphere.ncells - 1 # due to the mass stensil.
    ae        = grain_interface
    ac        = 0.5 * ae[1:] + 0.5 * ae[:-1]
    M         = sphere.M(sphere.rc)
    ## TODO sphere._vg[:,1] = (sphere._rc[:,2]-sphere._rc[:,1])/sphere.dt[1] APPLY THIS TO THE SPHERE GAS!
    fi             = func_dist(ae[:-1], ae[1:])/func_dist(ae[0], ae[-1])
    mi             = M[1:] - M[:-1]
    massdist       = mi[:,None] * fi[None,:] * D
    listsize       = len(fi) * len(mi)
    massdist.shape = listsize
    r              = np.repeat(sphere.rc[:-1], Ngrain)
    s              = np.tile(ac, Nbins)

    particles        = Particles(listsize, sphere, 0 ,1, 2,**kwards)
    particles.r[:,0] = r
    particles.v[:,:] = 0
    particles.s[:,0] = s
    particles.M      = massdist
    particles.D      = D
    particles.ngrains = Ngrain
    particles.sinterface = ae
    if debug:
       return particles, mi, fi
    return particles

class DataParticles(object):
    """ # rc, vg, rhog , M, T 
        #  0   1     2   3  4  """
        # NOTE only works on .npy files as of now (hardcoded)
        # TODO add time interpoaltion
        # NOTE TODO BUG i save my data when they are at the next step, ie at idstep = 0 i am actually at idstep = 1 in datasphere.
  #  __slots__ = '_rc', '_vg' ,'_rho', '_M' , '_T','_re', 'nsteps', 'ncells', 'idstep', 'dt', 't', 'interpolatevg',\
  #       'interpolaterho', 'interpolateM', 'interpolateT', 'mu','units' , 'interpolatesearch', '_dr'
    def __init__(self, loadpath, env, skiprows = 0, max_rows = None):
        self.dt = np.loadtxt(loadpath + '\\dt.txt'    , skiprows = skiprows,max_rows = max_rows)
        self.M  = np.loadtxt(loadpath + '\\mass.txt'  , skiprows = skiprows,max_rows= max_rows)
        self.r  = np.loadtxt(loadpath + '\\radius.txt', skiprows = skiprows,max_rows = max_rows)
        self.s  = np.loadtxt(loadpath + '\\s.txt'     , skiprows = skiprows,max_rows = max_rows)
        self.t  = np.loadtxt(loadpath + '\\time.txt'  , skiprows = skiprows,max_rows = max_rows)
        self.vg = np.loadtxt(loadpath + '\\vel.txt'   , skiprows = skiprows,max_rows = max_rows)
        self.env   = env
        self.units = env.units
        # interpolation arrays
    def search(self,rm):
        return np.int64(self.env.interpolatesearch[self.env.idstep](rm))
    
    def St(self,rm,s, snapid,rho_d = 1.6):
        
        r = np.zeros((len(rm),2))
        r[:,0] = rm
        rm     = magnitude(r,True)
        ts   = 1/self.env.drag(r,s,rho_d)
        tff  = np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm))) 
        return ts/tff

def dust_to_gas(p: DataParticles, s: float, snapid: int):
    """
    Returns the dust to gast ratio of the smallest particle.
    NOTE changes idstep

    input : (Union[Particles, DataParticles], int)
    -- Class of dataparticles or particles initialized with an envirment ; what snap id we should check
    TODO, remove snapid as an argument, then it worlds on other env end datasphere, given it has radial bins.

    return : (ndarray, ndarray, ndarray, ndarray)
    -- (dust to gas ratio, number of particles in bins, what bins where each particle is located, mask for smallest grain bin)
        Dust to gas format is D = [M_bin0, M_bin1, M_bin2...M_binM]
    """
    p.env.idstep = snapid
    s_low  = (p.s[snapid] == s)
    bins   = p.search(p.r[snapid,s_low])
    Mg     = p.env._M[:,snapid]
    dm_g   = Mg[1:] - Mg[:-1]
    dm_g   = dm_g
    
    m_dust_bin = [np.sum(p.M[snapid,s_low][bins == i]) for i in bins]
    N          = [np.sum(bins == i) for i in bins]
    binsrange  = np.arange(len(N))

    #if debug:
        #return 
    return m_dust_bin/dm_g[binsrange], N, bins, s_low, m_dust_bin, dm_g


class Particles_Corrected(Particles):
    """
    This here needs a new KDK that uses beta to correct the radius of the partilces radius in our interpolation functions, eg. if r_correct is found then we use self.r = self.r + dt*g(r_correct)
    TODO: ADD beta as a parameter for DataSphere as it only depends on the enviorment.
    TODO: IMPLEMENT and TEST inheritance from PARTICLES
    TODO: ADD KDK_drag_Corrected that uses the interpolation method that shows alpha and beta to calulate the corrected value.
    TODO: TEST that it works 

    TODO: Implement Grain distributed constructor for this class
    """
    def __init__(self, n, env, R, Stmin, Stmax,rho_d = 1.6, C = 0.004):
        super().__init__(n, env, R, Stmin, Stmax,rho_d, C)




    def r_trans(self):
        """
        a-based rtrans where f = (p-a)/(b-a)
        """
        idx = self.search
        rp  = magnitude(self.r)
    # if p < q
    #TODO NOTE due to how we do the velcotiy we dont have anything in the very last cell, and this can createa  problem due to dvg so we may need to limit it to it.
    #TODO espessially dvg
        rmask =   rp <= self.env.rc[idx] # TODO NOTE this changed this here to force beta to be 1 for rp < ra for idx == 0
        ## small hack for idx = ncells -1
        rup   = ( rp >  self.env.rc[idx]) & (idx == self.env.ncells-1)
        ###
        idx[(rmask) & (idx != 0)] -= 1# TODO may be wrong, small hack
        ## small hack for idx = ncells -1
        idx[rup] = -1
        ##
        #if idx == self.env.ncells-2:
            #raise('ErRRROROROR')
        # NOTE (idx != self.env.nsteps) this may not be needed in a version
        #      as we should always move down if we are located at the last cell
        Dt   = self.t[:,0] - self.env.t[self.env.idstep] # TODO assign this as a self. attribute
        Dr   = self.env.rc[1:] - self.env.rc[:-1] #TODO as above
        Dvg  = self.env._vg[1:,self.env.idstep] - self.env._vg[:-1,self.env.idstep] # NOTE here is a problem as the last cell ie dvg may be wrong
        
        
        beta = 1/(1+Dvg[idx]*Dt/Dr[idx]) 
        beta[(rmask) & (idx == 0)] = 1 # NOTE changed to one swe we just get a normal translation.
        ### small hack for idx = ncells -1
        beta[rup]                  = 1
        ####

        ra   = self.env.rc[idx]
        vga  = self.env._vg[idx, self.env.idstep]
    
        #res      = np.zeros_like(self.r)  
        res      = (beta*rp + (1-beta)*ra - Dt * beta * vga )[:,None] * self.r/magnitude(self.r, True)
        return res


    def kick(self,step=0.5):
        """
        Must be paired with self.v

        returns a step for v based on gravity        
        
        """
        force  = self.gravity()#not really a force, more so an acceleration
        self.v[self.runmask] = self.v[self.runmask] + step*self.dt[self.runmask]*force
        #return             step*self.dt*force
    
    def drift(self,step=1.0):
        """ 
        Must be paired with self.v

        returns a drift step for r based on velocity v.
        """
        self.r[self.runmask] = self.r[self.runmask] + step*self.dt[self.runmask]*self.v[self.runmask]
        #return             step*self.dt*self.v

    def gravity(self):
        rm = magnitude(self.r[self.runmask], True)
        return  - self.r[self.runmask] * self.units.G*self.env.M(rm)/(rm*rm*rm) 
        #return  - self.r * self.units.G*self.units.m_Sun/(rm*rm*rm) # TODO add mass change
        #return -self.r*self.G/(rm*rm*rm) # self .G is temp


    def KDK_drag_Corrected(self,same_dt=True):
        self.kick(0.5) # Update vel for grav
        rtrans            = self.r_trans()[self.runmask]
        r                 = rtrans  + 0.5 * self.dt[self.runmask]* self.v[self.runmask] # <----
        Cdt               = 0.5*self.dt[self.runmask] * self.env.drag(r,self.s[self.runmask],self.rho_d) 
        self.v[self.runmask] = self.v[self.runmask]  - Cdt*(self.v[self.runmask] - self.env.vg(r))/(1 + Cdt) 
        self.drift(1.0)
        self.kick(0.5)
        rtrans            = self.r_trans()[self.runmask] # this may be neeeded
        Cdt               = 0.5*self.dt[self.runmask] * self.env.drag(rtrans,self.s[self.runmask])
        self.v[self.runmask] = self.v[self.runmask]  - Cdt*(self.v[self.runmask] - self.env.vg(rtrans))/(1 + Cdt)

    def DKD_drag_Corrected(self,same_dt=True):
        #smask = (self.s == self.s.min())[:,0]

        self.drift(0.5) # Update vel for grav
        self.kick(1.0)
        rtrans = self.r_trans()
        r      = rtrans  + 0.5 * self.dt* self.v # <----
        Cdt    = 0.5*self.dt * self.env.drag(r,self.s,self.rho_d) 
        self.v = self.v  - Cdt*(self.v - self.env.vg(r))/(1 + Cdt) 
        self.drift(0.5)


    def KDK_drag_Corrected_collected(self,same_dt=True):
        r      = self.r + 0.5 * self.dt* self.v
        rm     = magnitude(r, True)
        dtgrav = 0.5 * self.dt * (- r * self.units.G*self.env.M(rm)/(rm*rm*rm))

        rtrans = self.r_trans()
        rt     = rtrans + 0.5 * self.dt* self.v # <----
        Cdt    = 0.5 * self.dt * self.env.drag(rt,self.s,self.rho_d) 
        self.v = self.v + ( dtgrav - Cdt*(self.v - self.env.vg(rt)) )/(1 + Cdt) 

        self.drift(1.0)

        rm     = magnitude(self.r, True)
        dtgrav   = 0.5*self.dt*(- self.r * self.units.G*self.env.M(rm)/(rm*rm*rm))

        rtrans = self.r_trans() # this may be neeeded
        Cdt    = 0.5*self.dt * self.env.drag(rtrans,self.s)
        self.v = self.v + ( dtgrav - Cdt*(self.v - self.env.vg(rtrans)) )/(1 + Cdt) 

    def courant_(self, Cdr = 0.2, Cff = 0.1, Csnap = 0.1, Ckep = 0.1, same_dt = False, verbose = False):# add verbose
        """
        Courant based on using the translatered radius to calulate dt.
        Seems to work better as it may more correclt correspond to where we actually pull our data.
        kinda like when we evolve.
        """
        r    = self.r_trans()
        dr   = self.env.dr[self.search,None]
        rm   = magnitude(r, True)
        a    = self.gravity()# - self.env.drag(r,self.s,self.rho_d)*(self.v - self.env.vg(r))# NOTE   g - c_d(vd - vg)#
        rr   = magnitude(self.r, True)
        if np.any(magnitude(a) == 0):
            #print('debug')
            self.gravity()# - self.env.drag(r,self.s,self.rho_d)*(self.v - self.env.vg(r))# NOTE   g - c_d(vd - vg)#


        dt_dr   = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        dt_ff   = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep  = Ckep  * np.sqrt(rm*rm*rm/(self.units.G*self.env.M(rm)))       # kepler based
        dt1     = np.minimum(dt_dr  , dt_ff) 
        dt2     = np.minimum(dt_snap, dt_kep)
        if np.any(dt_kep == np.inf):
            debug = True
        self.dt =  np.minimum(dt1, dt2)
        if verbose:
            imin = np.argmin( self.dt )
            print( 'idx: %i,  radius: %.3e, cell: %i, dt_dr: %.3e, dt_ff: %.3e, dt_snap: %.3e, dt_kep: %.3e, size: %.3e, a_grav: %.3e' \
                   %( imin, rr[imin]/self.units.AU, self.search[imin], dt_dr[imin]/self.units.yr,
                   dt_ff[imin]/self.units.yr, dt_snap/self.units.yr, dt_kep[imin]/self.units.yr , self.s[imin], a[imin][0]) )
        if same_dt:
            self.dt[:] = self.dt.min()

def grain_dist_radius_Correction(sphere:DataSphere, func_dist, grain_interface, D = 0.01, repeat = 1, debug = False, **kwards):
    """
    All starts along the positive x directions
    
    dust velocity is currently always 0 so same as gas to start.
    
    
    TODO This could be inplemented into the original such that we give the corrected class as an argument.


    args::
        sphere          : DataSphere
        -- The sphere which contrains interpolated data from a collapse
        func_dist       : callable
        -- Solution to the integral g(a1,a2) =  int_a1^a2 dn/dlna a^2 da
           used in fraction fj = g(a_l,a_u)/g(a_min,a_amx)
        grain_interface : 1D ndarray 
        -- array of the interfaces between grain bins 
        D               : float
        -- The dust-to-gas ratio, default = 1/100
           This defualt will give truncation error as 1/100 in binary is an infinite series.
    kwards::
        rho_d : float1
        -- material density of the dust, default = 1.6 g/cm^3
        C     : float
        -- Courant number, default = 0.004
    
    return: Particles Class
    -- The paricles class with grains created from a distribution
       The data is structure so [{a[:],radial bin 0},{a[:],radial bin 1},{a[:],radial bin 2},...] ie grains in radial bin 0 first then all grains in radial bin 2 etc.

    TODO BUG there seems to be a bug in placement of particles at the lower end, there are decimal diffrence 
    """
    # To add more grains of the same size i just need to repeat every stat, and divide the mass by the number of repeats.
    # ie i have [[s1,s2,s3], [s1,s2,s3]] --repeat 3 times--> [[s1,s1,s1,s2,s2,s2,s3,s3,s3], [s1,s1,s1,s2,s2,s2,s3,s3,s3]]
    # do this for all traits then do M = M/3 in this case, because we split each macro particels into 3
    Ngrain    = len(grain_interface) - 1 # number of grain bins not interfaces
    Nbins     = sphere.ncells - 1 # due to the mass stensil.
    ae        = grain_interface
    ac        = 0.5 * ae[1:] + 0.5 * ae[:-1]
    M         = sphere.M(sphere.re) # NOTE used to be sphere.rc
    ## TODO sphere._vg[:,1] = (sphere._rc[:,2]-sphere._rc[:,1])/sphere.dt[1] APPLY THIS TO THE SPHERE GAS!
    fi             = func_dist(ae[:-1], ae[1:])/func_dist(ae[0], ae[-1]) # fj from report
    mi             = M[1:] - M[:-1] # cell mass
    listsize       = len(fi) * len(mi)
    #massdist       = mi[:,None] * fi[None,:] * D
    #massdist.shape = listsize
    r              = np.repeat(sphere.re[:-1], Ngrain)#np.repeat(sphere.rc[:-1], Ngrain)
    dr             = np.repeat(sphere.dr, Ngrain) # should follow the evolution of r
    s              = np.tile(ac, Nbins)
    fj             = np.tile(fi, Nbins )
    ii             = np.tile(np.arange(1,repeat+1), listsize) # to work this has to start from 1

    ## code for mutliple grains of same size pr bin.
    listsize *= repeat
    dr        = np.repeat(dr      , repeat)
    r         = np.repeat(r       , repeat)
    r         = r + dr*(ii-0.5)/repeat # smoothly distributes the particles
    s         = np.repeat(s       , repeat)
    fj        = np.repeat(fj, repeat)


    ## correcting mass with rho scaling
    rho      = sphere.rho(r)#<---
    dr_part  = dr/repeat
    dV       = 4/3*np.pi * ( (r + 0.5*dr_part)**3 - (r - 0.5*dr_part)**3 )
    massdist = D * fj * rho * dV
    #massdist  = np.repeat(massdist, repeat)/repeat
    Ngrain  *= repeat
    ## \---
    particles        = Particles_Corrected(listsize, sphere, 0 ,1, 2,**kwards)
    particles.r[:,0] = r
    particles.v[:,:] = 0
    particles.s[:,0] = s
    particles.M      = massdist
    particles.D      = D # dust to gas ratio number.
    particles.ngrains = Ngrain
    particles.sinterface = ae # TODO make this a get setter that gives the grain bins
    return particles, fj, D

def dust_to_gas_ParticlesClass(p: DataParticles,size : float):
    """
    Returns the dust to gast ratio of the smallest particle.
    NOTE changes idstep

    input : (Union[Particles, DataParticles], int)
    -- Class of dataparticles or particles initialized with an envirment ; what snap id we should check
    TODO, remove snapid as an argument, then it worlds on other env end datasphere, given it has radial bins.

    return : (ndarray, ndarray, ndarray, ndarray)
    -- (dust to gas ratio, number of particles in bins, what bins where each particle is located, mask for smallest grain bin)
        Dust to gas format is D = [M_bin0, M_bin1, M_bin2...M_binM]
    """
    s_low  = (p.s == size)[:,0]
    bins   = p.search[s_low]
    Mg     = p.env._M[:,p.env.idstep]
    dm_g   = Mg[1:] - Mg[:-1]
    
    #m_dust_bin = [np.sum(p.M[s_low][bins == i]) for i in bins] #
    m_dust_bin = np.zeros_like(dm_g)
    farray     = np.zeros_like(dm_g)
    for pm, r, sbin in zip(p.M[s_low], magnitude(p.r)[s_low], bins):
        #if sbin == 813 or sbin-1 == 813 or sbin +1 == 813:
        #    print('deb')
        #if sbin == len(dm_g)-1:
        #    print('deb')
        if r <= p.env.rc[sbin]:
            if sbin == 0:
                m_dust_bin[sbin] = m_dust_bin[sbin] + pm
            else:
                dr = p.env.rc[sbin]-p.env.rc[sbin-1] #center to center
                f  = ( p.env.rc[sbin] - r )/dr
                m_dust_bin[sbin]     = m_dust_bin[sbin] + (1-f)*pm
                m_dust_bin[sbin - 1] = m_dust_bin[sbin - 1] + f*pm
        else:
            if sbin == p.env.ncells - 1:
                m_dust_bin[sbin] = m_dust_bin[sbin] + pm
            else:
                dr = p.env.rc[sbin + 1]-p.env.rc[sbin] #center to center
                f  = ( r - p.env.rc[sbin] )/dr
                m_dust_bin[sbin]     = m_dust_bin[sbin] + (1-f)*pm
                m_dust_bin[sbin + 1] = m_dust_bin[sbin + 1] + f*pm


    N          = [np.sum(bins == i) for i in bins]
    binsrange  = np.arange(len(N))

    #if debug:
        #return 
    # m_dust_bin/dm_g[binsrange]
    return m_dust_bin/dm_g, N, bins, s_low, m_dust_bin, dm_g



def dust_to_gas_DataParticles(p,size : float, snapid: int):
    """
    Returns the dust to gast ratio of the smallest particle.
    NOTE changes idstep

    input : (Union[Particles, DataParticles], int)
    -- Class of dataparticles or particles initialized with an envirment ; what snap id we should check
    TODO, remove snapid as an argument, then it worlds on other env end datasphere, given it has radial bins.

    return : (ndarray, ndarray, ndarray, ndarray)
    -- (dust to gas ratio, number of particles in bins, what bins where each particle is located, mask for smallest grain bin)
        Dust to gas format is D = [M_bin0, M_bin1, M_bin2...M_binM]
    """
    sphereid_old = p.env.idstep
    debugid = 0
    p.env.idstep = snapid
    s_low  = (p.s[snapid] == size)
    rm     = p.r[snapid]
    bins   = p.search(rm)[s_low]
    Mg     = p.env._M[:,snapid]
    dm_g   = Mg[1:] - Mg[:-1]
    rc     = p.env.re[:-1] + 0.5*p.env.dr  # p.env._rc # rc     = p.env._re[:-1,:] + 0.5*p.env.dr[:,None] 
    #m_dust_bin = [np.sum(p.M[s_low][bins == i]) for i in bins] #
    m_dust_bin = np.zeros_like(dm_g)
    farray     = np.zeros_like(dm_g)
    Mlost      = 0
    for pm, r, sbin in zip(p.M[snapid, s_low], p.r[snapid, s_low], bins):
        if (sbin == p.env.ncells-1):
            Mlost += pm
            continue

        if r <= rc[sbin]:
            if sbin == 0:
                #dr = rc[sbin, snapid]-rc[sbin - 1,snapid] #center to center
                #f  = ( rc[sbin, snapid] - r )/dr
                #m_dust_bin[sbin]     = m_dust_bin[sbin]     + (1-f)*pm
                m_dust_bin[sbin] = m_dust_bin[sbin] + pm
                f = 0
            else:
                dr = rc[sbin]-rc[sbin - 1] #center to center
                f  = ( rc[sbin] - r )/dr
                m_dust_bin[sbin]     = m_dust_bin[sbin]     + (1-f)*pm
                m_dust_bin[sbin - 1] = m_dust_bin[sbin - 1] + f*pm
                ####
                farray[sbin-1] = farray[sbin-1] + f
                ####
        else:
            if sbin == p.env.ncells - 2: # its -2 because p.env.ncells i counted, so it effectivly start from 1, so we have ncells-1
                                         # cells indexed from 0. as we reduce this by 1 cell in dmg, then we should check sbin p.en.ncells -2
                                         # eg. dm_g[4095] is out of bound but dm_g[4094] is the last index
                m_dust_bin[sbin] = m_dust_bin[sbin] + pm
                f = 0
            else:
                dr = rc[sbin + 1]-rc[sbin] #center to center
                f  = ( r - rc[sbin] )/dr
                m_dust_bin[sbin]     = m_dust_bin[sbin] + (1-f)*pm
                m_dust_bin[sbin + 1] = m_dust_bin[sbin + 1] + f*pm
                ####
                farray[sbin+1] = farray[sbin+1] + f
                ####
        farray[sbin] = farray[sbin] + (1-f)

    N            = [np.sum(bins == i) for i in bins]
    binsrange    = np.arange(len(N))
    p.env.idstep = sphereid_old
    #if debug:
        #return 
    # m_dust_bin/dm_g[binsrange]
    #return m_dust_bin/dm_g, N, bins, s_low, m_dust_bin, dm_g, farray
    return m_dust_bin/dm_g, m_dust_bin ,farray, rc, Mlost


def CalcDustmass_snapid(p : DataParticles,s,snapid):
    old_snap     = p.env.istep
    p.env.idstep = snapid
    smask        = (p.s[snapid] == s)
    f            = np.zeros(p.n, dtype = np.float64)
    neighbor     = np.zeros(p.n, dtype = np.int64  ) 
    dustmass_bin = np.zeros(p.n, dtype = np.float64)
    rc           = p.env._re[:-1,:] + 0.5*p.env.dr[:,None]  # p.env._rc
    rc           = rc[:,snapid]
    bins         = p.search(p.r)
    for idx, pm, r, sbin in enumerate(zip(p.M[snapid, smask], p.r[snapid, smask], bins)):
        if (r <= rc[sbin]):
            if (sbin == 0):
                f[idx]             = 0
                dustmass_bin[sbin] += pm
                neighbor[idx]      = sbin
            else:
                dr                      = rc[sbin] - rc[bin - 1]
                f[idx]                 = (rc[sbin] - r)/dr
                dustmass_bin[sbin]     += (1 - f[idx])*pm 
                dustmass_bin[sbin - 1] += f[idx]*pm
                neighbor[idx]          = sbin -1
        else:
            if (sbin == p.ncells - 2):
                f[idx]             = 0
                dustmass_bin[sbin] += pm
                neighbor[idx]      = sbin
            else:
                dr                      = rc[sbin + 1] - rc[sbin]
                f[idx]                 = (r - rc[sbin])/dr
                dustmass_bin[sbin]     += (1 - f[idx])*pm 
                dustmass_bin[sbin - 1] -= f[idx]*pm
                neighbor[idx]          = sbin + 1
    p.env.idsnap = snapid_old
    return dustmass_bin, f, neighbor


def CalcDustmass_size(p : DataParticles, s):
    smask        = (p.s == s).flatten()
    f            = np.zeros(p.n           , dtype = np.float64)
    neighbor     = np.zeros(p.n           , dtype = np.int64  ) 
    dustmass_bin = np.zeros(p.env.ncells-1, dtype = np.float64) # should be cell- 2 long
    rc           = p.env.re[:-1] + 0.5*p.env.dr # p.env._rc
    bins         = p.search
    rm           = magnitude(p.r)
    # TODO NOTE f[smask][idx] is temp, i should be f[idx] when i seperate dustmass_bin from this func
    F_temp       = np.zeros_like(f[smask])
    neigh_temp   = np.zeros_like(neighbor[smask])
    for idx, (pm, r, sbin) in enumerate(zip(p.M[smask], rm[smask], bins[smask])):
        if (r <= rc[sbin]):
            if (sbin == 0):
                F_temp[idx]         = 0
                dustmass_bin[sbin] += pm
                neigh_temp[idx]     = sbin
            else:
                dr                      = rc[sbin] - rc[sbin - 1]
                F_temp[idx]             = (rc[sbin] - r)/dr
                dustmass_bin[sbin]     += (1 - F_temp[idx])*pm 
                dustmass_bin[sbin - 1] += F_temp[idx]*pm
                neigh_temp[idx]         = sbin -1
        else:
            if (sbin == p.env.ncells - 2):
                F_temp[idx]         = 0
                dustmass_bin[sbin] += pm
                neigh_temp[idx]     = sbin
            else:
                dr                      = rc[sbin + 1] - rc[sbin]
                F_temp[idx]             = (r - rc[sbin])/dr
                dustmass_bin[sbin]     += (1 - F_temp[idx])*pm 
                dustmass_bin[sbin + 1] += F_temp[idx]*pm
                neigh_temp[idx]        = sbin + 1
    f[smask]        = F_temp 
    neighbor[smask] = neigh_temp
    return dustmass_bin, f, bins, neighbor  # NOTE F is temp


def CalcDustmass(p : DataParticles):#, s):
    f            = np.zeros(p.n           , dtype = np.float64)
    neighbor     = np.zeros(p.n           , dtype = np.int64  ) 
    #dustmass_bin = np.zeros(p.env.ncells-1, dtype = np.float64) # should be cell- 2 long
    rc           = p.env.re[:-1] + 0.5*p.env.dr # p.env._rc
    bins         = p.search
    rm           = magnitude(p.r)
    # TODO NOTE f[smask][idx] is temp, i should be f[idx] when i seperate dustmass_bin from this func
    #F_temp       = np.zeros_like(f[smask])
    #neigh_temp   = np.zeros_like(neighbor[smask])
    for idx, (pm, r, sbin) in enumerate(zip(p.M, rm, bins)):
        if (r <= rc[sbin]):
            if (sbin == 0):
                f[idx]           = 0
                #dustmass_bin[sbin] += pm
                neighbor[idx]    = sbin
            else:
                dr                      = rc[sbin] - rc[sbin - 1]
                f[idx]                  = (rc[sbin] - r)/dr
                #dustmass_bin[sbin]     += (1 - F_temp[idx])*pm 
                #dustmass_bin[sbin - 1] += F_temp[idx]*pm
                neighbor[idx]           = sbin -1
        else:
            if (sbin == p.env.ncells - 2):
                f[idx]              = 0
                #dustmass_bin[sbin] += pm
                neighbor[idx]       = sbin
            else:
                dr                      = rc[sbin + 1] - rc[sbin]
                f[idx]                  = (r - rc[sbin])/dr
                #dustmass_bin[sbin]     += (1 - F_temp[idx])*pm 
                #dustmass_bin[sbin + 1] += F_temp[idx]*pm
                neighbor[idx]           = sbin + 1
    #f[smask]        = F_temp 
    #neighbor[smask] = neigh_temp
    return f, bins, neighbor  # NOTE F is temp


def InitialMassCorrection(p, fj, D, max_iter=100 ,diff_eps = 1e-15, verbose = False):
    """
    fj  is the rho fraction for each grian size,
    D is the dust to gas ratio. default is 0.01 ie 1/100
    TODO start with max_iter, ie set number of iterations
    TODO start with a given size, then do all sizes
    TODO f, bins, neighbor are all constant for a grain speciecs (maybe for all grains).
         This means that ones calulated, they dont need to be updated, so we can update cell_mass indepnedly?
    """   
    Mg        = p.env._M[:, p.env.idstep]
    gas_mass = Mg[1:] - Mg[:-1]
    ## cell_mass should be repeated so each particles get the correct mass 
    # to iterate over
    #f1, bins1, neighbor1 = CalcDustmass(p)

    s_unique = np.unique(p.s)
    cell_mass1 = np.zeros( len(gas_mass) )
    if verbose:
        print('Correcting mass of %i of unique grains'%len(s_unique))
        Mold = p.M.sum()
    for ids, ss in enumerate(s_unique):
        cell_mass, f, bins, neighbor = CalcDustmass_size(p, ss)
        d2g = cell_mass/gas_mass
        k = 0 
        while  k < max_iter: # ((d2g.max() - d2g.min()) > diff_eps) 
            cell_mass1[:] = 0
            smask = (p.s == ss).flatten()

            FF_bin       = D*fj*gas_mass[bins]     / cell_mass[bins]
            FF_neigh     = D*fj*gas_mass[neighbor] / cell_mass[neighbor]
            p.M[smask]   = ( ( 1-f[smask] )*FF_bin[smask] + f[smask]*FF_neigh[smask] )*p.M[smask]
            cell_mass, f, bins, neighbor = CalcDustmass_size(p, ss)
            d2g = cell_mass/gas_mass

            k += 1
            #### TODO Attepmt at a more optimized version, doesn't work 
            #own   = (1-f1[smask])*p.M[smask]
            #neigh = (f1[smask])*p.M[smask]
            #np.add.at(cell_mass1, bins1[smask], own)
            #np.add.at(cell_mass1, neighbor1[smask], neigh)
            ####


            #bmask = bins    [smask]
            #nmask = neighbor[smask
            #FF_bin1       = D*fj[bmask]*gas_mass[bmask] / cell_mass1[bmask]
            #FF_neigh1     = D*fj[nmask]*gas_mass[nmask] / cell_mass1[nmask] 

            #p.M[smask]   = ( ( 1-f1[smask] )*FF_bin1 + f1[smask]*FF_neigh1 )*p.M[smask]
        if verbose:
            print('S = %.2e'%ss,' cm  ','Iter:',k,' d2g array:', d2g)
    if verbose:
        print('Old total mass: %.4e' %(Mold/p.units.m_Sun), "m_sun ,New total mass: %.4e"%( p.M.sum()/p.units.m_Sun),' m_sun')
        print('Done correcting dust mass')

class save_particles(object):
    def __init__(self, p: Particles, savepath: str):
        """
        needs exit enter metods,then i should be albe to use with

        If n as input this with quit the entire script!!
        """
        self.particles = p
        self.savepath  = savepath
        try: os.mkdir(self.savepath)
        except FileExistsError:
            inputt = input('File Folder exists, files might be interfered with, want to continure (y/n)').lower()
            if inputt != 'y':
                print('Closing script')
                quit()#

        self.files  = ['dt.txt', 'mass.txt', 'radius.txt',
                       's.txt',  'time.txt',    'vel.txt',
                       'x.txt',     'y.txt',    'vx.txt', 
                       'vy.txt',   'St.txt', 'rtrans.txt',
                       'dragx.txt', 'gravityx.txt']
        self.open()

    def open(self):
        """
        Opens and potentially creates creates the folders
        """
        self.filehandle = {}
        for savefile in self.files:
            self.filehandle[savefile] = open(self.savepath + savefile, 'w')
            ## instead of append we might need write? then we wont add to existing data
            ## just overwrite it
      
    def save(self):
        """
        Append saves to all
        """
        p = self.particles
        np.savetxt(self.filehandle['dt.txt']      , p.dt.T          , newline = ' ')
        np.savetxt(self.filehandle['mass.txt']    , p.M             , newline = ' ') 
        np.savetxt(self.filehandle['s.txt']       , p.s             , newline = ' ') 
        np.savetxt(self.filehandle['time.txt']    , p.t             , newline = ' ') 
        np.savetxt(self.filehandle['radius.txt']  , magnitude(p.r)  , newline = ' ') 
        np.savetxt(self.filehandle['x.txt']       , p.r[:,0]        , newline = ' ') 
        np.savetxt(self.filehandle['y.txt']       , p.r[:,1]        , newline = ' ')
        np.savetxt(self.filehandle['vel.txt']     , magnitude(p.v)  , newline = ' ') 
        np.savetxt(self.filehandle['vx.txt']      , p.v[:,0]        , newline = ' ')
        np.savetxt(self.filehandle['vy.txt']      , p.v[:,1]        , newline = ' ')
        np.savetxt(self.filehandle['St.txt']      , p.St            , newline = ' ')
        np.savetxt(self.filehandle['rtrans.txt']  , p.r_trans()[:,0], newline = ' ')
        gravity = p.gravity()[:,0] 
        np.savetxt(self.filehandle['gravityx.txt'], gravity         , newline = ' ')
        


        # writes "enter" or newline into each file
        for handle in self.filehandle.values():
            handle.write('\n')    

    def close(self):
        """
        Closes all 
        should it delete self?
        """
        for savefile in self.filehandle.values():
            savefile.close()
        del self


def main():
    ##-------##
    kep_percent = 0
    savefolder = './SimData/'
    foldername = 'Both_Grains_4_pr_5bins_1000_%.2f_kep_2/'%(kep_percent)

    ##-------##

    sphere  = DataSphere('./run_extracted/run067/',sc.CGS)  
    # OLD ? sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
    
    sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
    sphere.update()

    ginterface  = np.geomspace(1.42921456e-06, 1.42921456e-01, 6)
    #ginterface  = np.geomspace(1.42921456e-06, 1.42921456e-04, 6)
    p,fj,D       = grain_dist_radius_Correction(sphere, GrainDistIntegrated, ginterface,debug = True, repeat = 4)
    p.v[:,1]    = kep_percent*np.sqrt(p.units.G*p.env.M(magnitude(p.r))/magnitude(p.r)) # NOTE assumes that p.v is ineritrely 0
    InitialMassCorrection(p, fj, D,max_iter = 100,verbose = True)
    filehandler = save_particles(p, savefolder+foldername) 
    filehandler.save()
    print('Beginning Simulation')
    for idt, dtstep in enumerate(p.env.dt[:-1],1):
        k=0
        if idt ==40:
            debug = True

        # set mask true thus i can use it in while.
        while np.any(p.runmask):
            p.courant(Csnap = 0.001,Ckep = 0.1,verbose = False)
            p.dt = np.minimum(p.env.t[p.env.idstep + 1] - p.t, p.dt)
            p.runmask[(p.t >= p.env.t[idt]).flatten()] = False
            #p.[p.t >= p.env.t[idt]] = 0
            p.KDK_drag_Corrected()
            p.t = p.t + p.dt # <--- NOTE should be after kdk
            k = k + 1
            if k == 1000000:
                raise('Max iter it, returning')
        p.env.next
        p.runmask[:] = True
        filehandler.save()
        print("snap %i:   \titer = %i"%(idt,k))
        ##if idt ==180:
        ##    filehandler.close()
        ##    quit()
        #if idt == 1:
        #    filehandler.close()
        #    quit() 
    print('Done at t = %.2f yr'%(np.mean(p.t)/p.units.yr))
    filehandler.close()

def microns_to_cm(x):
    return x/10000

def sim(folder, kep, spherefolder):
    ##-------##
    kep_percent = kep
    ns          = 6
    s0, s1      = microns_to_cm(21), microns_to_cm(117)  # 1.42921456e-03, 1.42921456e-01 # <- already in cm
    si,sj       = microns_to_cm(25), microns_to_cm(100)
    vkep1,vkep2 = 0.025, 1  # kep

    savefolder = './SimData/'
    #foldername = folder + kep
    foldername = 'Both_Grains_4_pr_5bins_1000_%.2fkep/'%(kep)
    print(foldername)
    
    ##-------##
    vkep_dist = scilate.interp1d([si,sj],[vkep1, vkep2], kind = 'linear', bounds_error=False, fill_value=(vkep1, vkep2))
    ##-------##
    sphere  = DataSphere('./run_extracted/' + spherefolder +'/',sc.CGS)  
    # OLD ? sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
    
    sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
    sphere.update()

    ginterface  = np.geomspace(s0, s1, ns)
    #ginterface  = np.geomspace(1.42921456e-06, 1.42921456e-04, 6)
    p,fj,D      = grain_dist_radius_Correction(sphere, GrainDistIntegrated, ginterface,debug = True, repeat = 4)
    ##### for a dist of vkep. uncommoent for constant
    #kep_percent = vkep_dist(p.s).flatten()
    #####
    p.v[:,1]    = kep_percent*np.sqrt(p.units.G*p.env.M(magnitude(p.r))/magnitude(p.r))
    InitialMassCorrection(p, fj, D,max_iter = 100,verbose = True)
    filehandler = save_particles(p, savefolder+foldername + '/') 
    filehandler.save()
    print('Beginning Simulation')
    for idt, dtstep in enumerate(p.env.dt[:-1],1):
        k=0
        if idt ==3:
            debug = True

        # set mask true thus i can use it in while.
        while np.any(p.runmask):
            p.courant(Csnap = 0.001,Ckep = 0.1,verbose = False)
            p.dt = np.minimum(p.env.t[p.env.idstep + 1] - p.t, p.dt)
            p.runmask[(p.t >= p.env.t[idt]).flatten()] = False
            #p.[p.t >= p.env.t[idt]] = 0
            p.KDK_drag_Corrected()
            p.t = p.t + p.dt # <--- NOTE should be after kdk
            k = k + 1
            if k == 1000000:
                raise('Max iter it, returning')
        p.env.next
        p.runmask[:] = True
        filehandler.save()
        print("snap %i:   \titer = %i"%(idt,k))
        ##if idt ==180:
        ##    filehandler.close()
        ##    quit()
        #if idt == 1:
            #filehandler.close()
        #    quit() 
        
    print('Done at t = %.2f yr'%(np.mean(p.t)/p.units.yr))
    filehandler.close()

if __name__ == '__main__':
    #sim('run086_25-100_large_',0.9,"run086")
    #sim('run074_25-100_microns_','Dist',"run074")

    #sim('run067_large_Grains_4_pr_5bins_1000_', 0.9, "run086")
   # sim('run067_large_Grains_4_pr_5bins_1000_', 0.7, "run086")
    sim('run067_large_Grains_4_pr_5bins_1000_', 0.5, "run086")
    #sim('run067_large_Grains_4_pr_5bins_1000_', 0.3, "run086")
    #sim('run067_large_Grains_4_pr_5bins_1000_', 0.1, "run086")
    #sim('run067_large_Grains_4_pr_5bins_1000_', 0.0, "run086")



    #sim('large_Grains_4_pr_5bins_1000_', 0.9)





