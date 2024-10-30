import os
import numpy as np

import scipy.integrate as scint
import scipy.optimize as sciopt
import scipy.interpolate as scilate


def numerical_generation(xstart = 1e-6, xend = 8, n = 64, write_file = False):
    """
    Calulates the critical values and the adimmenional radius and density.
    """
    def lane_emden(y, xi):
        """
        This is not the adimensional vairable from the text!
        y = [u, z]
        dydt = [dudx,dzdx]
        """
        u, z = y
        dydx = [z/(xi*xi), -xi*xi*np.exp(u)]
        return dydx
    xi, dx = np.linspace(xstart, xend, n, retstep= True)
    y_init = [0, 0]
    res    = scint.odeint(lane_emden, [0, 0], xi)
    yi     = np.exp(res[:,0])
    m      = np.cumsum(xi*xi*yi*dx)
    # critical values regeneration #
    mc      = np.sqrt(yi)*m
    dmdx    = np.gradient(mc,xi)
    dmdy    = np.gradient(mc,yi)
    mc_f    = scilate.interp1d(xi,mc,'linear')
    dmdx_f  = scilate.interp1d(xi,dmdx,'linear')
    dmdy_f  = scilate.interp1d(yi,dmdy,'linear')
    xcrit   = sciopt.fsolve(dmdx_f,6)
    ycrit   = sciopt.fsolve(dmdy_f,1/14.1)#1/1.41 seems to give an error?
    mcrit   = mc_f(xcrit)
    Acrit   = mcrit/(np.sqrt(ycrit)*xcrit)
    crit    = np.hstack([Acrit, xcrit, ycrit, mcrit])
    if (write_file == True):
        try:
            os.mkdir('./emden-solve')
        except FileExistsError:
            pass
        np.savetxt('./emden-solve/crit-value.txt', crit)
        np.savetxt('./emden-solve/rho-adim.txt'  , np.c_[xi,yi,m])
    else:
        return crit,np.c_[xi,yi,m].T

def extract_data(path: str, savetype: str = 'npy', extractpath:str = './run_extracted/'):
    """
     # rc, vg, rhog , M, T 
     #  0   1     2   3  4
     this is for [:,:]

    savetype : npy or txt
    """
    # TODO add a save folder argument
    #    t  dt
    # 0  0  10 <-- this dt is to go from 0 to 1 etc 
    # 1 10  9  
    #    .  .
    #    .  .
    # n tn  0  <--- this dt should be 0 as we are done 
    if path[-1] != '/': path = path + '/'
    savepath    = extractpath + path[2:]
    try:     os.stat (extractpath)
    except:  os.mkdir(extractpath)
    try:     os.stat (savepath)
    except:  os.mkdir(savepath)

    output_files = os.listdir(path)

    data   = np.loadtxt(path + output_files[0],usecols = (0,2,6,8,10,1)) # NOTE first file handled here for pre assigning memory
    nsteps = len(output_files)
    ncells = len(data[:,0])
    time   = np.zeros((2,nsteps))
    with open(path + output_files[0]) as filedata: # NOTE t for would always start at 0 so this is just in case
        tstring     = filedata.readline()
        ieq         = tstring.find('=') + 1 ; isec = tstring.find('s')
        time[0,0]   = np.float64(tstring[ieq:isec]) 
    data_extracted  = np.zeros((ncells, 6, nsteps))
    data_extracted[:,:, 0] = data
    for idx, output in enumerate(output_files[1:],1):
        with open(path + output) as filedata:
            tstring     = filedata.readline()
            ieq         = tstring.find('=') + 1 ; isec = tstring.find('s')
            time[0,idx] = np.float64(tstring[ieq:isec]) 
            data_extracted[:,:,idx] = np.loadtxt(filedata,usecols = (0,2,6,8,10,1))
    time[1,:-1] = time[0, 1:] - time[0, :-1]
    if  savetype == 'npy':
        np.save(savepath + 'data',data_extracted)
        np.save(savepath + 'time', time)
    elif savetype == 'txt':
        np.savetxt(savepath + 'data',data_extracted)
        np.savetxt(savepath + 'time', time)
    print('Done Extrating')
#extract_data('./run067')
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
    for pm, r, sbin in zip(p.M[s_low], magnitude(p.r)[s_low], bins):
        if r <= p.env.rc[sbin]:
            if sbin == 0:
                m_dust_bin[sbin] = m_dust_bin[sbin] + pm
            else:
                dr = p.env.rc[sbin]-p.env.rc[sbin-1] #center to center
                f  = ( p.env.rc[sbin] - r )/dr
                m_dust_bin[sbin]     = m_dust_bin[sbin] + (1-f)*pm
                m_dust_bin[sbin - 1] = m_dust_bin[sbin - 1] + f*pm
        else:
            if sbin == p.env.nsteps - 1:
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
    return m_dust_bin/dm_g[binsrange], N, bins, s_low, m_dust_bin, dm_g


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
    for pm, r, sbin in zip(p.M[snapid][s_low], p.r[snapid][s_low], bins):
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

