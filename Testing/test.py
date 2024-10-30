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

class Particles(object):
    __slots__ = 'n','t', 'dt' , 'env','units','r','s', 'Stmin' , 'Stmax', 'rho_d', 'v','M' ,'D', 'ngrains','sinterface'
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

    def courant_velocity(self,C = 0.01,same_dt=True):
        """ Determine time step(s) for at least self.nstep steps per orbit.
            If same_dt is True, take the same time step for all particles,
            If False, take nstep steps per orbit, for all radii
        """
        vm             = magnitude(self.v)
        gm             = magnitude(self.gravity())
        self.dt[:,0]   = C*vm/gm
        #self.dt[np.isnan(self.dt)] = 0 # TODO NOTE remember this is here.
        if same_dt:
            self.dt[:] = self.dt.min()
        # self.t = self.t + self.dt

    def courant(self, Cdr = 0.2, Cff = 0.1, Csnap = 0.1, Ckep = 0.1, same_dt = False):
        dr   = self.env.dr[self.search,None]
        rm   = magnitude(self.r, True)
        a    = self.gravity() - self.env.drag(self.r,self.s,self.rho_d)*(self.v - self.env.vg(self.r))# NOTE   g - c_d(vd - vg)

        #dt_dr   = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        dt_ff   = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep  = Ckep  * np.sqrt(rm*rm*rm/(self.units.G*self.env.M(rm)))       # kepler based
        dt1     = dt_ff#np.minimum(dt_dr, dt_ff) 
        dt2     = np.minimum(dt_snap,dt_kep)
        self.dt =  np.minimum(dt1, dt2)
        if same_dt:
            self.dt[:] = self.dt.min()
        

    def gravity(self):
        rm = magnitude(self.r, True)
        return  - self.r * self.units.G*self.env.M(rm)/(rm*rm*rm) # TODO add mass change
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
        self.kick(0.5)
        # NOTE + != * :(
        r      = self.r  # + 0.5 * self.dt* self.v # half here as it needs to sync with kick
        Cdt    = 0.5*self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
        self.v = self.v  - Cdt*(self.v - self.env.vg( r ))/(1 + Cdt) # vg = 0
        self.drift(1.0)
        self.kick(0.5)
        Cdt    = 0.5*self.dt * self.env.drag(self.r,self.s)
        self.v = self.v  - Cdt*(self.v - self.env.vg(self.r))/(1 + Cdt) # vg = 0



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
         'interpolaterho', 'interpolateM', 'interpolateT', 'mu','units' , 'interpolatesearch', '_dr'
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
        self._dr     = data[1:,5,:] - data[:-1,5,:] 
        self.idstep  = idstep # what sna
        self.dt      = time[1]
        self.t       = time[0]
        self.mu      = mu
        self.units   = units
        # interpolation arrays
        self.update()

    def update(self):
        n = range(self.ncells)
        self.interpolatevg     = [scilate.interp1d(self._rc[:,idx], self._vg[:,idx] ,fill_value = (self._vg [0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolaterho    = [scilate.interp1d(self._rc[:,idx], self._rho[:,idx],fill_value = (self._rho[0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolateM      = [scilate.interp1d(self._rc[:,idx], self._M[:,idx]  ,fill_value = (self._M  [0,idx],self._M[-1,idx]), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolateT      = [scilate.interp1d(self._rc[:,idx], self._T[:,idx]  ,fill_value = (self._T  [0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolatesearch = [scilate.interp1d(self._re[:,idx], n               ,fill_value = (0               ,self.ncells - 1), bounds_error = False) for idx in range(self.nsteps)]
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
         return self.interpolatevg[self.idstep](rmag)

    def vg(self, r):
        """
        TODO NOTE this needs to be a 2d vector, in cartesian corodiantes, this is the radial
        """
        rm = magnitude(r,True)
        vm = self.interpolatevg[self.idstep](rm) # NOTE this is signed
        return vm * r/rm

    def rho(self, rmag):
        return self.interpolaterho[self.idstep](rmag)
    
    def M(self, rmag):
        return self.interpolateM[self.idstep](rmag)
    
    def T(self, rmag):
        return self.interpolateT[self.idstep](rmag)

    # Copied from sphere class
    def drag(self, r, s,rho_d = 1.6):
        """
        This is actually the inverse stopping time.
        so drag = 1/t_s
        rho_d = 1.6 g/cm^3
        """
        rm    = magnitude(r,True)
        rho_g = self.rho(rm)
        vth   = np.sqrt(8/np.pi * self.units.k_B*self.T(rm)/(self.mu*self.units.m_p))
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


    def r_trans_b(self):
        """
        b-based rtrans where f = (b-p)/(b-a)
        TODO NOTE proerbly use this as we should have any particles over sphere.rc[-1]
                  but we do have particles that will get smaller that sphere.rc[0]
        """
        idx = self.search
        rp  = magnitude(self.r)

        idx[rp < self.env.rc[idx]] -= 1

        Dt   = self.t[:,0] - self.env.t[self.env.idstep]
        Dr   = self.env.rc[1:] - self.env.rc[:-1]
        Dvg  = self.env._vg[1:, self.env.idstep] - self.env._vg[:-1, self.env.idstep]

        beta = 1/(1 + Dvg[idx]*Dt/Dr[idx])##
        
        rb  = self.env.rc[idx+1]
        vgb = self.env._vg[idx+1, self.env.idstep]

        res      = np.zeros_like(self.r)
        res[:,0] = beta * rp + (1 - beta)*rb - beta* Dt * vgb
        return res


    def r_trans(self):
        """
        a-based rtrans where f = (p-a)/(b-a)
        """
        idx = self.search
        rp  = magnitude(self.r)
    # if p < q
    #TODO NOTE due to how we do the velcotiy we dont have anything in the very last cell, and this can createa  problem due to dvg so we may need to limit it to it.
    #TODO espessially dvg
        idx[(rp < self.env.rc[idx]) & (idx != 0)] -= 1# TODO may be wrong, small hack
        #if idx == self.env.ncells-2:
            #raise('ErRRROROROR')
        # NOTE (idx != self.env.nsteps) this may not be needed in a version
        #      as we should always move down if we are located at the last cell
        Dt   = self.t[:,0] - self.env.t[self.env.idstep] # TODO assign this as a self. attribute
        Dr   = self.env.rc[1:] - self.env.rc[:-1] #TODO as above
        Dvg  = self.env._vg[1:,self.env.idstep] - self.env._vg[:-1,self.env.idstep] # NOTE here is a problem as the last cell ie dvg may be wrong
        beta = 1/(1+Dvg[idx]*Dt/Dr[idx]) 

        ra   = self.env.rc[idx]
        vga  = self.env._vg[idx, self.env.idstep]
        
        res      = np.zeros_like(self.r)
        res[:,0] = beta*rp + (1-beta)*ra - Dt * beta * vga   #temp hack, needs to be this but needs a r/|r| so it scales correctly.
        return res
        

    def KDK_drag_Corrected(self,same_dt=True):
        smask = (self.s == self.s.min())[:,0]
        rtrans = self.r_trans()
        self.kick(0.5) # Update vel for grav
        r      = rtrans  + 0.5 * self.dt* self.v # <----
        Cdt    = 0.5*self.dt * self.env.drag(r,self.s,self.rho_d) 
        self.v = self.v  - Cdt*(self.v - self.env.vg(r))/(1 + Cdt) 
        self.drift(1.0)
        self.kick(0.5)
        rtrans = self.r_trans() # this may be neeeded
        Cdt    = 0.5*self.dt * self.env.drag(rtrans,self.s)
        self.v = self.v  - Cdt*(self.v - self.env.vg(rtrans))/(1 + Cdt) 


    def courant_(self, Cdr = 0.2, Cff = 0.1, Csnap = 0.1, Ckep = 0.1, same_dt = False, verbose = False):# add verbose
        """
        Courant based on using the translatered radius to calulate dt.
        Seems to work better as it may more correclt correspond to where we actually pull our data.
        kinda like when we evolve.
        """
        r    = self.r_trans()
        dr   = self.env.dr[self.search,None]
        rm   = magnitude(r, True)
        a    = self.gravity() - self.env.drag(r,self.s,self.rho_d)*(self.v - self.env.vg(r))# NOTE   g - c_d(vd - vg)#
        rr   = magnitude(self.r, True)


        dt_dr   = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        dt_ff   = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep  = Ckep  * np.sqrt(rr*rr*rr/(self.units.G*self.env.M(rm)))       # kepler based
        dt1     = np.minimum(dt_dr  , dt_ff) 
        dt2     = np.minimum(dt_snap, dt_kep)
        self.dt =  np.minimum(dt1, dt2)
        if verbose:
            imin = np.argmin( self.dt )
            print( 'idx: %i,  radius: %.3e, cell: %i, dt_dr: %.3e, dt_ff: %.3e, dt_snap: %.3e, dt_kep: %.3e, size: %.3e, a_grav: %.3e' \
                   %( imin, rr[imin]/self.units.AU, self.search[imin], dt_dr[imin]/self.units.yr,
                   dt_ff[imin]/self.units.yr, dt_snap/self.units.yr, dt_kep[imin]/self.units.yr , self.s[imin], a[imin][0]) )
        if same_dt:
            self.dt[:] = self.dt.min()

def grain_dist_radius_Correction(sphere:DataSphere, func_dist, grain_interface,D = 0.01,debug = False, **kwards):
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

    particles        = Particles_Corrected(listsize, sphere, 0 ,1, 2,**kwards)
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



def dust_to_gas_DataParticles(p: DataParticles,size : float, snapid: int):
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
    s_low  = (p.s[snapid] == size)
    rm     = p.r[snapid]
    bins   = p.search(rm)[s_low]
    Mg     = p.env._M[:,p.env.idstep]
    dm_g   = Mg[1:] - Mg[:-1]
    
    #m_dust_bin = [np.sum(p.M[s_low][bins == i]) for i in bins] #
    m_dust_bin = np.zeros_like(dm_g)
    farray     = np.zeros_like(dm_g)
    for pm, r, sbin in zip(p.M[snapid,s_low], magnitude(p.r)[snapid,s_low], bins):
        #if sbin == 813 or sbin-1 == 813 or sbin +1 == 813:
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


def main():
    folder = 'Datarun_dragdt'
    path   = './SimData/'+folder
    try:  os.mkdir(path)
    except FileExistsError:
        answer = input('File exist, do you want to continue (y/n)')
        if answer.lower() != 'y':
            raise('Folder Exists')
    fdt     = open(path + '\\dt.txt','a')
    fmass   = open(path + '\\mass.txt','a')
    fradius = open(path + '\\radius.txt','a')
    fs      = open(path + '\\s.txt','a')
    ftime   = open(path + '\\time.txt','a')
    fvel    = open(path + '\\vel.txt','a')
    sphere  = DataSphere('./run_extracted/run067/',sc.CGS)  
    
    sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
    sphere.update()

    ginterface = np.geomspace(1.42921456e-06, 1.42921456e-01, 20)

    p,_,_ = grain_dist_radius_Correction(sphere, GrainDistIntegrated, ginterface,debug = True)
    #dtg, _, _ ,_ ,_,_= dust_to_gas_ParticlesClass(p,p.s.min())
    print('RUN!')
    np.savetxt(fdt    , p.dt.T        , newline=' '); fdt.write('\n')
    np.savetxt(fmass  , p.M           , newline=' '); fmass.write('\n')
    np.savetxt(fradius, magnitude(p.r), newline=' '); fradius.write('\n')
    np.savetxt(fs     , p.s           , newline=' '); fs.write('\n')
    np.savetxt(ftime  , p.t           , newline=' '); ftime.write('\n')
    np.savetxt(fvel   , magnitude(p.v), newline=' '); fvel.write('\n')
    for idt, dtstep in enumerate(p.env.dt[:-1],1):
        k = 0

        while np.any(p.t < p.env.t[idt]):
            
            p.courant_(Csnap = 0.1,verbose = False)
            p.dt = np.minimum(p.env.t[idt] - p.t, p.dt)
            p.dt[p.t >= p.env.t[idt]] = 0
            p.KDK_drag_Corrected()
            p.t = p.t + p.dt # <--- NOTE should be after kdk
            k = k + 1
            if k == 1000000:
                raise('Max iter it, returning')

        p.env.next
        print('snap: %i done \t iterations here: %i'%(idt, k)) 
        np.savetxt(fdt    , p.dt.T        , newline=' '); fdt.write('\n')
        np.savetxt(fmass  , p.M           , newline=' '); fmass.write('\n')
        np.savetxt(fradius, magnitude(p.r), newline=' '); fradius.write('\n')
        np.savetxt(fs     , p.s           , newline=' '); fs.write('\n')
        np.savetxt(ftime  , p.t           , newline=' '); ftime.write('\n')
        np.savetxt(fvel   , magnitude(p.v), newline=' '); fvel.write('\n')
    #if idt == 1:
        #    dtg, _, _ ,_ ,_,_= dust_to_gas_ParticlesClass(p,p.s.min())
        print("snap %i:\titer = %i"%(idt,k))
    print('Done at t = %.2f yr'%(np.mean(p.t)/p.units.yr))
    fdt.close()    
    fmass.close()  
    fradius.close()
    fs.close()     
    ftime.close()  
    fvel.close()  

if __name__ == '__main__':
    main()