import numpy as np

import scipy.integrate as scint
import scipy.optimize as sciopt
import scipy.interpolate as scilate

import scaling as sc
from .ultils import magnitude

class BE_sphere(object):
    # xmin, xmax are temp? maybe NOTE
    __slots__ = 'units','rho_c','rc','xcrit','y', 'xmin','xmax','mu','rho','T','M','rmin', 'R','vg'
               #sc.cgs , float, float, float, callable, float, float
    def __init__(sphere,units, M_BE = 1.0, T = lambda x: 10,drag=0.001, vg = lambda x: 0):
        Acrit,xcrit,ycrit,mcrit = np.loadtxt('./emden-solve/crit-value.txt',unpack = True) 
        xi,yi,_                 = np.loadtxt('./emden-solve/rho-adim.txt'  ,unpack = True)
        dx = xi[1]-xi[0]

        M_BE  = M_BE * units.m_Sun
        cs_sq = units.k_B*T/(units.mu*units.m_p) # base units
        rho_c = mcrit**2 * cs_sq**3/(4*np.pi*ycrit * M_BE**2 * units.G**3)
        r_c   = np.sqrt(cs_sq/(4*np.pi * units.G * rho_c))
        if callable(T): sphere.T = T
        else:           sphere.T = lambda x : T
        
        sphere.units = units
        sphere.rc    = r_c
        sphere.rho_c = rho_c
        sphere.mu    = units.mu


        sphere.xcrit = xcrit
        sphere.xmin  = xi[0]    ; sphere.xmax = xi[-1]
        sphere.rmin  = xi[0]*r_c; sphere.R    = xcrit*r_c

        m = 4*np.pi*np.cumsum(xi*xi * yi)*dx*rho_c * r_c**3
        sphere.M     = scilate.interp1d(xi*r_c,m)
        sphere.y     = scilate.interp1d(xi,yi,'linear')
        sphere.rho   = scilate.interp1d(xi*r_c,yi*rho_c,'linear')

        r  = np.linspace(sphere.rmin, sphere.R,len(xi))
        m  = sphere.M(r) ; rho = sphere.rho(r)
        ## To ensure that that are no out of bounds errors, this menas that M is constant out of R, rho is 0 as the sphere is stable and is equal to R
        sphere.M   = scilate.interp1d(r, m, 'linear',  fill_value = (m[0],m[-1]), bounds_error = False)
        sphere.rho = scilate.interp1d(r, rho,'linear', fill_value = (rho[0], 0), bounds_error = False) # unsure here if it should be zero
        sphere.vg  = vg
        if callable(vg): sphere.vg = vg
        else:            sphere.vg = lambda x : vg
        #fill_value = (m[0],m[-1]), bounds_error = False
        #fill_value = (yi[0],yi[-1]), bounds_error = False
        #fill_value= (yi[0]*rho_c,yi[-1]*rho_c), bounds_error = False

    def drag(self, r, s,rho_d = 1.6,):
        """
        This is actually the inverse stopping time.
        so drag = 1/t_s


        rho_d = 1.6 g/cm^3
        """
        rm    = magnitude(r,True)
        rho_g = self.rho(rm)
        vth   = np.sqrt(8/np.pi * self.units.k_B*self.T(rm)/(self.mu*self.units.m_p))
        return rho_g*vth/(rho_d*s)
        #vth =  self.units.k_B*self.T/(self.units.mu*self.units.m_p)

    # def drag_coeff(self,r,particle_size):
        # rm = magnitude(r)
        # return self.drag*self.density(rm)/particle_size

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

        dt_dr   = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        dt_ff   = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep  = Ckep  * np.sqrt(rm*rm*rm/(self.units.G*self.env.M(rm)))       # kepler based
        dt1     = np.minimum(dt_dr, dt_ff) 
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
        idx[rp < self.env.rc[idx]] -= 1

        Dt   = self.t[:,0] - self.env.t[self.env.idstep] # TODO assign this as a self. attribute
        Dr   = self.env.rc[1:] - self.env.rc[:-1] #TODO as above
        Dvg  = self.env._vg[1:,self.env.idstep] - self.env._vg[:-1,self.env.idstep]
        beta = 1/(1+Dvg[idx]*Dt/Dr[idx]) 

        ra   = self.env.rc[idx]
        vga  = self.env._vg[idx, self.env.idstep]
        
        res      = np.zeros_like(self.r)
        res[:,0] = beta*rp + (1-beta)*ra - Dt * beta * vga   #temp hack, needs to be this but needs a r/|r| so it scales correctly.
        return res
        

    def KDK_drag_Corrected(self,same_dt=True):
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


    def courant_(self, Cdr = 0.2, Cff = 0.1, Csnap = 0.1, Ckep = 0.1, same_dt = False):
        """
        Courant based on using the translatered radius to calulate dt.
        Seems to work better as it may more correclt correspond to where we actually pull our data.
        kinda like when we evolve.
        """
        r    = self.r_trans()
        dr   = self.env.dr[self.search,None]
        rm   = magnitude(r, True)
        a    = self.gravity() - self.env.drag(r,self.s,self.rho_d)*(self.v - self.env.vg(r))# NOTE   g - c_d(vd - vg)#

        dt_dr   = Cdr   * np.sqrt(2*dr/magnitude(a,True))                       # NOTE a is the total acceleration ie dv/dt that is solve for velocity
        dt_ff   = Cff   * np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rm)))   # free fall time.
        dt_snap = Csnap * self.env.dtsnap                                       # Snapshot based 
        dt_kep  = Ckep  * np.sqrt(rm*rm*rm/(self.units.G*self.env.M(rm)))       # kepler based
        dt1     = np.minimum(dt_dr, dt_ff) 
        dt2     = np.minimum(dt_snap,dt_kep)
        self.dt =  np.minimum(dt1, dt2)
        if same_dt:
            self.dt[:] = self.dt.min()









