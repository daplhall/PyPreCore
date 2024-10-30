# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# ## Todo
# Fra Troels:
# Her er noget som ville være spændende at få sat op:
# 
# 1) Lav et simpelt setup for en prestellar core. Brug en Bonnor-Ebert sphere (se evt comp astro uge 4 eller uge 5)
# 
# 2) Antag at BE-sfæren er i ligevægt, og derfor ikke ændrer sige
# 
# 3) Beregn hvordan støv vil udvikle sig, når det slippes løs i sådan i sfære, som funktion af støv størrelsen.
# 
# Brug gerne rigtige tal. Dvs. R_BE ≈  10,000 AU, M_BE ≈2 Msun, T = 10 K
# 
# 
# Mig selv:
# 
# 1) For inskrevet at man kan gøre det fysisk. Brug de critiske ligninger osv. du har fundet
# 
# 2) Test om det passer med massen.
# 
# 3) implementer partilker. udnyt at det minder om comp astro examen, (start med et spestemt støv regime)
# 
# 4) implementer drag. 
# 
# 5) test at ting ikke flyver af helvede til.
# 
# 6) find på en måde at teste hvor gyldigt det er. (stopping time?)
# 
# 7) Ryd op i kode.
# 
# bonus:
#     implementer scaling pakken
# %% [markdown]
# ## Excat solution to Bonnar-emden sphere through numerical analsys.
# (Used Rober Estalella, The Bonnor-Ebert Sphere, provied in the computational astophyscis course)
# 
# If we assume that we have a homogenous gas, in all directions except the radial. This gas i also assumed to be in hydrostatic equlibrium.
# \begin{equation}
#     \frac{dP}{dr} = - \frac{G M_r}{r^2}\rho
# \end{equation}
# Here $M_r$ is the mass, which depends on our radius r and can be discrubed by the density.
# \begin{equation}
#     dM_r = 4\pi r^2 \rho dr \Leftrightarrow \frac{dM_r}{dr} = 4\pi r^2\rho
# \end{equation}
# We can now isolate M_r i our hydrostatic equlibrium equation and insert it into our mass equation.
# \begin{equation}
#     \frac{1}{r^2}\left(\frac{r^2}{\rho}\frac{dP}{dr}\right) + 4\pi G\rho = 0
# \end{equation}
# We can express the preasure in terms of the sound speed and the density by assuming a isothermal gas ie. $P = c_s^2 \rho$, where $c_s^2 = k_b T/(\mu m_h)$.
# \begin{equation}
#     \frac{1}{r^2}\left(\frac{r^2}{\rho}\frac{d\rho}{dr}\right) + \frac{4\pi G}{c_s^2}\rho = 0
# \end{equation}
# The boundry conditions for this system is 
# \begin{align}
#     \rho(0) &= \rho_c\\
#     \frac{d\rho}{dr}|_{r=0} &= 0
# \end{align}
# This is because we want the center of the sphere to be physcially representitive with a finite center density, and we don't want it to in crease in negative r, as this is not physcial.
# 
# Before we solve the system it would be beneficial to rewrite it with to a dimmensionalless form. Here we need to define some charateristik parameters. In this equation we only have units of density and length, so a good choice of charateristik parameters would be $\rho_c$ which we define ourself, and $r_c$ which we need to find though dimmensional analasys. $r_c$ can be found though dimmensional analsys.
# \begin{equation}
#     r_c = \frac{c_s}{\sqrt{4\pi G \rho_c}}
# \end{equation}
# $4\pi$ is inclueded as it is dimmensionalless and it is beneficial to include it, specific for our equation, as it allows us to get a nice adimmensional equation.
# 
# We can now define adimmensional parameters x,y.
# \begin{align}
#     x = \frac{r}{r_c}\\
#     y = \frac{\rho}{\rho_c}
# \end{align}
# By inserting these into our density equation we get a adimmensial equation of form.
# \begin{equation}
#     \frac{1}{x^2} \frac{d}{dx}\left(\frac{x^2}{y}\frac{dy}{dx}\right) + y = 0
# \end{equation}
# We can make this reseble the Lane–Emden equation by defining $u = \ln(y)$ and iserting it into our adimmensional equation. This is done as we get a every easy equation to solve numerically.
# \begin{equation}
#     \frac{1}{x^2}\frac{d}{dx}\left(x^2 \frac{du}{dx}\right) + e^u = 0
# \end{equation}
# This ODE can be writte as a set of coulpled ODEs, that makes it a lot easier to solve, by defining $z = x^2 d_xu$
# \begin{equation}
#     \frac{du}{dx} = \frac{z}{x^2}
# \end{equation}
# \begin{equation}
#     \frac{dz}{dx} = - x^2 e^u
# \end{equation}
# Now we just need to define the boundries in these new units.
# \begin{align}
#     u(0) &= \ln(y(0)) = \ln\left(\frac{\rho_c}{\rho_c}\right) = 0\\
#     z(0) &= 0^2\frac{du}{dx}|_0 = 0
# \end{align}
# Now we can solve this system numerially with a forward method like runge-kutta.
# %% [markdown]
# ## The Mass
# The adimensional bonnar-emden for our boundry conditions, is somewhat independent (we need to solve further or equal to our bonnar sphere) from our physcial representation. So to be able to set up a system of a BE-sphere of a given a mass and radius we need to couple some of our charateristic parameters with the mass $M_r$.
# 
# If we start with our radial mass equation
# \begin{equation}
#     dM_r = 4\pi r^2 \rho dr
# \end{equation}
# If make this adimensional and then integrate both sides from the centor the the edge then we might be able to get an idea of the central density as the mass should be total at $R_{max}$.
# \begin{align}
#     dM_r &= 4\pi r_c^3 \rho_c \; x^2 y dx\\
#     M_r  &=  4\pi r_c^3 \rho_c \int_0^x x'^{\;2} y dx'
# \end{align}
# This integral equation can be estimated numerically by writing it as a sum, if we assume that we have a small mesh spacing. One could also use the trapzoidal rule, but as we define our own mesh we just use the sum discretization.
# \begin{equation}
#     M_r = 4\pi r_c^3 \rho_c \sum_{i=0}^n x_i^2y_i \Delta x_i = 4\pi r_c^3 \rho_c m_c
# \end{equation}
# where $x_n = x$ ie. the upper bound of our integral. $m_c$ is a adimensional mass, that is depended on the adimensional laden-emden equation.
# 
# Important note: $r_c$ depends on $\rho_c$
# So write
# \begin{equation}
# M_r = 4 \pi \left(\frac{c_s}{\sqrt{4\pi G \rho_c}}\right)^3 \rho_c m_c \Leftrightarrow
# \end{equation}
# %% [markdown]
# ## Critical stability values
# To have a stable cloud we need
# \begin{equation}
#     \frac{\partial P_0}{\partial r_0} = c_s^2 \frac{\partial \rho_0}{\partial r_0} < 0
# \end{equation}
# For if this was > 0 for the outer region then the gas would diffuse out. So we can try and find the criticla mass value for a given chouse of $\rho(r)$ in this case $\rho_0$. remember $y = \rho/\rho_c$
# \begin{equation}
#  M_r(r) = \frac{c_s^3}{\sqrt{4\pi \rho_c} G^{3/2} } m_c(x) = \frac{c_s^3\sqrt{y}}{\sqrt{4\pi \rho(r)} G^{3/2} } m_c(x) = \frac{c_s^3}{\sqrt{4\pi \rho} G^{3/2} } m'_c(x)
# \end{equation}
# We can now use $m'_c$ to find where the critical point for a cloud is, by discretizing it the same way as $m_c$. Here we need to find the first maxima, as this is this will give us the minimum radius that the cloud is allowed for it to be stable, this we can do in both the adimmensional space of x and y. When these are found we can find the mass and radius at these points simply by converting our critical adimmensional parameters.
# \begin{align}
#     r(r_o)= R_{BE} &= x_{crit} r_c = x_{crit} \frac{c_s}{\sqrt{4\pi G \rho_c}}\\
#     M(r_o) = M_{BE} &= m_{crit} \frac{c_s^3}{\sqrt{4\pi \rho_o} G^{3/2} } = m_{crit} \frac{c_s^3}{\sqrt{4\pi y_{crit} \rho_c} G^{3/2} }
# \end{align} 
# A relation between $M_{BE}, R_{BE} and T$ can be found by writing $M_{BE}$ in terms of $c_s$ and $R_{BE}$
# \begin{equation}
#     M_{BE} = \frac{m_{crit}}{\sqrt{y_{crit}}x_{crit}} \frac{c_s^2}{G} R_{BE} = A_{crit}\frac{c_s^2}{G} R_{BE}
# \end{equation}
# This we can write in practical units and we can let $\mu$ and $A_{crit}$ be a tunable parameter as it is unitless. $A_{crit}$ is a constant, but i will load it in from a file for sanity, plus i will allow recalulation of solving the adimmensional density, again as a sanity check.
# \begin{equation} 
#     \left[\frac{M_{BE}\;R_{BE}^{-1}}{M_\odot \; AU^{-1}}\right] = 9.303038506e\text{-}6\frac{A_{crit}}{\mu}\left[\frac{T}{K}\right] \; \text{where} \; A_{crit} \approx 2.43
# \end{equation}
# We can now define our charatistic values like the charatistic length $r_c$ and the center denisty $\rho_c$ from a given combination of 2 of mass, radius and temperture.
# \begin{equation}
# M_{BE} = \frac{m_{crit}}{\sqrt{y_{crit}}} \frac{c_s^3}{\sqrt{4\pi \rho_c} G^{3/2} } \Leftrightarrow  \rho_c= \frac{m_{crit}^2}{y_{crit}} \frac{c_s^6}{4\pi M_{BE}^2 G^{3} }
# \end{equation}
# %% [markdown]
# ## Particle gas drag
# To find the drag of a particle we need to look at aerodynamic drag forces (Philip J. Armitage lecture notes).
# 
# \begin{equation}
#     F_D = -\frac{1}{2}C_d\pi s^2\rho_g v_{dg}^2
# \end{equation}
# for the Epstein regine we have that
# \begin{equation}
#     C_D = \frac{8}{3}\frac{\bar{v}}{v_{dg}} \mspace{10mu} \text{for} \mspace{10mu} s < \frac{9}{4}\lambda_f
# \end{equation}
# Here $\lambda_f$ is the mean free path, and $\bar{v} = \sqrt{8/\pi}\mspace{5mu}c_s$ (<- bar v).
# This gives us a drag force of 
# \begin{equation}
#     F_D = -\frac{V_d}{s}\rho_g \bar{v} v_{dg}
# \end{equation}
# Here $V_d$ is the volume of a given dust grain.
# Here we can simply write the acceleration from the drag and such form newtons second law.
# \begin{equation}
#     m_d a = F_{total} = F_d+\sum_i F_i \Leftrightarrow a = \frac{F_d}{m_d} + \frac{1}{m_d}\sum_i F_i = -\frac{V_d}{m_d}\frac{\bar{v}}{s}\rho_g v_{dg}+\frac{1}{m_d}\sum_i F_i
# \end{equation}
# So here we get
# \begin{equation}
# a = -\frac{\rho_g}{\rho_d}\frac{\bar{v}}{s}v_{dg}+\frac{1}{m_d}\sum_i F_i
# \end{equation}
# TODO WRITE THE UPDATE STENSIL/SCHEME, implicit!
# 
# 
# The stopping time for this drag can easily be found through dimmensional analysis. Here we take the momentum change of the drag-dust velocity and the drag force.
# \begin{equation}
#     t_{stop} = \frac{|p_{dg}|}{|F_D|} = \frac{m_d v_{dg}}{\frac{V_d}{s}\rho \bar{v} v_{dg}} = \frac{\rho_d}{\rho_g}\frac{s}{\bar{v}}
# \end{equation}
# This is concequently the inverse of the method \lstinline{drag} for \lstinline{class Particles}.
# 
# From the stopping time and free-fall time
# \begin{equation}
#     t_{ff} = \sqrt{\frac{3\pi}{32G\rho_g}}
# \end{equation}
# 
# The stokes number is defined as 
# \begin{equation}
#     S_t = \frac{t_{stop}}{t_{ff}} = \frac{\rho_d s}{\rho_g(r) \bar{v}}\sqrt{\frac{32 G \rho_g}{3\pi}}
# \end{equation}
# So as $c_s$ is constant due to the sphere being isothermal, we have that the only r dependcy of $S_t$ is $\rho_g$.
# \begin{equation}
#     S_t \propto \frac{1}{\rho_g(r)}
# \end{equation}
# This means that in the inner parts of the sphere, our stokes number only depends on the grain size.
# To convert from stokes number to given size at a given radius we use.
# \begin{equation}
#     s = \frac{t_{ff} \rho_g \bar{v}}{\rho_d} St
# \end{equation}
# %% [markdown]
# ## Sphere based on snapshots.
# 
# To use af snapshot approach to the spheres evolution one has to take into account a few things. The first is the integration method.
# Here each snapshot has a given time given to them named  $t_n = [ t_0, t_1, t_2, ..., t_N]$ where $N$ is the number of snapshots. So we need to integrate our particles up till we hit the snapshot time, ie. start at $t_0$ then evolve particles, until the evolved time $t_evo == t_1$ and so forth.
# To control the step size one just choses the minimum between the particales courant conditions or the remaining time ie.
# \begin{equation}
#     \Delta t = \text{min}(\Delta t_{courant}, t_n - t_{evo})
# \end{equation}
# 
# This is repeated until all partiles are at $t_n$.
# 
# ASCI illustration
# \begin{enumerate}
# 1  •------>        | next iteration   •------>----->  | iterate until all  •------>----->->| Change to new snapshot
# 2  •-------------->|  same snapshot   •---------------| is at the next     •---------------| and repeat procces
# 3  •---------->    | ------------>    •---------->--->| shot.              •---------->--->| until no snapshots
# 4  •- >            |                  •->--->         |                    •->--->--->---->|   
# \end{enumerate}

# %%
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

try:
    os.stat( './Plots')
except: 
    os.mkdir('./Plots')
    os.mkdir('./Plots/drag')


# TODO make these into rc params.
labeloptions = {
    'fontsize' : 14     ,
    'fontname' : 'serif',
}
plotoptions = {
    'linewidth' : 2,
}

legendsoptions = {
    #'loc'      : 2      ,
    'frameon'  : False  ,

    'prop'     : {'family': 'serif',
                  'size'  : 12},
}


# %%
def magnitude(r, keepdims = False):
    """ 
    scalar magnitude of a vector
    from num astro
    """
    return np.sqrt(np.sum(r**2,1,keepdims=keepdims))

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
    @profile
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
    


# %%
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


# %%
class Particles(object):
    __slots__ = 'n','t', 'dt' , 'env','units','r','s', 'Stmin' , 'Stmax', 'rho_d', 'C', 'v'
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
        self.C      = C
        st     = np.geomspace((Stmin,), (Stmax,), n)
        rmag   = magnitude(self.r, True)
        tff    = np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rmag)))                           # TODO move to enviorment
        vth    = np.sqrt(8/np.pi * self.units.k_B*self.env.T(rmag)/(self.env.mu*self.units.m_p)) # TODO move to enviorment instead, and then just call it.
        self.s = tff * self.env.rho(rmag)* vth / self.rho_d * st

    @property
    def St(self):
        rmag = magnitude(self.r, True)
        ts   = 1/self.env.drag(self.r,self.s,self.rho_d)
        tff  = np.sqrt(3*np.pi/(32*self.units.G*self.env.rho(rmag))) 
        return ts/tff
    @profile
    def courant_velocity(self,same_dt=True):
        """ Determine time step(s) for at least self.nstep steps per orbit.
            If same_dt is True, take the same time step for all particles,
            If False, take nstep steps per orbit, for all radii
        """
        vm             = magnitude(self.v)
        gm             = magnitude(self.gravity())
        self.dt[:,0]   = self.C*vm/gm
        if same_dt:
            self.dt[:] = self.dt.min()
        # self.t = self.t + self.dt

    def courant_size(self,same_dt=True):
        """ Determine time step(s) for at least self.nstep steps per orbit.
            If same_dt is True, take the same time step for all particles,
            If False, take nstep steps per orbit, for all radii
        """
        s            = p.s[:,0]
        gm           = magnitude(self.gravity())
        self.dt[:,0] = self.C * np.sqrt(s/gm)
        if same_dt:
            self.dt[:] = self.dt.min()
        # self.t = self.t + self.dt
    @profile
    def gravity(self):
        rm = magnitude(self.r, True)
        M  = self.env.M(rm)
        return  - self.r * self.units.G*M/(rm*rm*rm) # TODO add mass change
        #return  - self.r * self.units.G*self.units.m_Sun/(rm*rm*rm) # TODO add mass change
        #return -self.r*self.G/(rm*rm*rm) # self .G is temp
    @profile
    def kick(self,step=0.5):
        """grav velocity step"""
        force  = self.gravity()#not really a force, more so an acceleration
        self.v = self.v + step*self.dt*force
    
    def drift(self,step=1.0):
        """ Spacial step"""
        self.r = self.r + step*self.dt*self.v

    def KDK(self,same_dt=True):
        self.courant_velocity()
        self.kick(0.5)
        self.drift(1.0)
        self.kick(0.5)
    @profile
    def KDK_drag(self,same_dt=True):
        self.kick(0.5)
        # NOTE + != * :(
        r      = self.r  + 0.5 * self.dt* self.v # half here as it needs to sync with kick
        Cdt    = self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
        self.v = self.v  - Cdt*(self.v - self.env.vg(r))/(1 + Cdt) # vg = 0
        self.drift(1.0)
        self.kick(0.5)
        Cdt    = self.dt * self.env.drag(self.r,self.s)
        self.v = self.v  - Cdt*(self.v - self.env.vg(self.r))/(1 + Cdt) # vg = 0


# %%
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

    data   = np.loadtxt(path + output_files[0],usecols = (0,2,6,8,10)) # NOTE first file handled here for pre assigning memory
    nsteps = len(output_files)
    ncells = len(data[:,0])
    time   = np.zeros((2,nsteps))
    with open(path + output_files[0]) as filedata: # NOTE t for would always start at 0 so this is just in case
        tstring     = filedata.readline()
        ieq         = tstring.find('=') + 1 ; isec = tstring.find('s')
        time[0,0]   = np.float64(tstring[ieq:isec]) 
    data_extracted  = np.zeros((ncells, 5, nsteps))
    data_extracted[:,:, 0] = data

    for idx, output in enumerate(output_files[1:],1):
        with open(path + output) as filedata:
            tstring     = filedata.readline()


            
            ieq         = tstring.find('=') + 1 ; isec = tstring.find('s')
            time[0,idx] = np.float64(tstring[ieq:isec]) 
            data_extracted[:,:,idx] = np.loadtxt(filedata,usecols = (0,2,6,8,10))
    time[1,:-1] = time[0, 1:] - time[0, :-1]
    if  savetype == 'npy':
        np.save(savepath + 'data',data_extracted)
        np.save(savepath + 'time', time)
    elif savetype == 'txt':
        np.savetxt(savepath + 'data',data_extracted)
        np.savetxt(savepath + 'time', time)
#extract_data('./run067')


# %%
class DataSphere(object):
    """ # rc, vg, rhog , M, T 
        #  0   1     2   3  4  """
        # NOTE only works on .npy files as of now (hardcoded)
        # TODO make the parameters here interpolations, ie interpolate over cells for a given time.
        # TODO add function that changes snap shot. Is this needed?
    __slots__ = '_rc', '_vg' ,'_rho', '_M' , '_T', 'nsteps', 'ncells', 'idstep', 'dt', 't', 'interpolatevg', 'interpolaterho', 'interpolateM', 'interpolateT', 'mu','units'
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
        self.idstep  = idstep
        self.dt      = time[1]
        self.t       = time[0]
        self.mu      = mu
        self.units   = units
        # interpolation arrays
        self.interpolatevg  = [scilate.interp1d(self._rc[:,idx], self._vg[:,idx] ,fill_value = (self._vg [0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolaterho = [scilate.interp1d(self._rc[:,idx], self._rho[:,idx],fill_value = (self._rho[0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolateM   = [scilate.interp1d(self._rc[:,idx], self._M[:,idx]  ,fill_value = (self._M  [0,idx],self._M[-1,idx]), bounds_error = False) for idx in range(self.nsteps)]
        self.interpolateT   = [scilate.interp1d(self._rc[:,idx], self._T[:,idx]  ,fill_value = (self._T  [0,idx],              0), bounds_error = False) for idx in range(self.nsteps)]

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

# %% [markdown]
# # Plotting
# ### Non collapsing BEsphere

# %%
def main():
    fig_rt  = plt.figure(figsize=(16,19))
    fig_rr  = plt.figure(figsize=(16,19))
    angles  = np.array([30])
    M, N    = 3,2
    sphere  = BE_sphere(sc.CGS,1,10)

    runtime  = 112817*sphere.units.yr #np.sqrt(3*np.pi/(32*sphere.units.G*sphere.rho(sphere.R))) 
    print('runtime:  t_ff(r_outer) =  %i yr'% (runtime/sphere.units.yr))
    for idx, angle in enumerate(angles):
        p        = Particles(10,sphere,R= sphere.R,Stmin= 1e-4,Stmax = 1e1)
        p.C      = 0.1
        pr       = 1 # %
        v_Kepler = np.sqrt(p.units.G*p.env.M(p.env.R)/p.env.R) * pr
        p.v[:,0] = -v_Kepler*np.sin(angle *np.pi/180)
        p.v[:,1] = +v_Kepler*np.cos(angle *np.pi/180)
        
        St = (p.St)#.astype(str)
        r  = []
        t  = []

        r.append(p.r)
        t.append(p.t)
        while np.all(p.t < runtime):
            p.courant_velocity(same_dt = True)
            p.dt[p.t >= runtime] = 0
            p.t = p.t + p.dt
            p.KDK_drag()
            r.append(p.r)
            t.append(p.t)

        r  = np.array(r); t = np.array(t)
        rr = np.sqrt(np.sum(r**2,2))
        leg = ['%.2e'%(q) for q in St]
        ax = fig_rt.add_subplot(M, N ,idx+1)
        ax.set_title('angle = %i'%(angle))
        ax.plot(t[:,:,0]/runtime,rr/p.units.AU)
        ax.set_xlabel(r'time $[t_{ff} = %.2e$ $yr]$'%(runtime/p.units.yr))                
        ax.set_ylabel(r'radius [AU]')
        ax.legend(leg,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_rt.tight_layout()
        # ------- r-r plots ------
        ax2 = fig_rr.add_subplot(M, N ,idx+1)
        ax2.set_title('angle = %i'%(angle))
        ax2.plot(r[:,:,0]/p.units.AU,r[:,:,1]/p.units.AU)
        ax2.set_xlabel(r'rad')                
        ax2.set_ylabel(r'radius [AU]')
        ax2.set_xlim([-p.env.R/p.units.AU ,p.env.R/p.units.AU])
        ax2.set_ylim([-p.env.R/p.units.AU ,p.env.R/p.units.AU])
        ax2.legend(leg,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_rr.tight_layout()
        #------------------------
        print('angle %i done' %(angle))
    fig_rt.savefig('./Plots/drag/%.2f-r-t.png'%(pr))
    fig_rr.savefig('./Plots/drag/%.2f-r-r.png'%(pr))  

if __name__ == '__main__':
    main()


# %%


# %% [markdown]
# ### Using data from A grid of 1D low-mass star formation collapse models, N. Vaytet & T. Haugbølle 

# %%
def main():
    global sphere
    fig_rt  = plt.figure(figsize=(16,19))
    fig_rr  = plt.figure(figsize=(16,19))
    fig_Stt = plt.figure(figsize=(16,19))

    angles  = np.array([0, 15, 30, 45, 60, 75])
    M, N    = 3,2
    sphere  = DataSphere('./run_extracted/run067/',sc.CGS)
    R       = sphere.rc[-1]
    runtime = sphere.t[-1]
    rct = np.repeat(sphere._rc[-1,:],2)[:-1]
    trc = np.repeat(sphere.t, 2)[1:]
    print('runtime:  t_{collapse} =  %i yr'% (runtime/sphere.units.yr))
    for idx, angle in enumerate(angles): # nan kommer af ændring i angle
        sphere.reset
        p        = Particles(10,sphere,R= R,Stmin= 1e-4,Stmax = 1e1)
        p.C      = 0.1
        pr       = 1 # %
        v_Kepler = np.sqrt(p.units.G*p.env.M(R)/R) * pr
        p.v[:,0] = -v_Kepler*np.sin(angle *np.pi/180)
        p.v[:,1] = +v_Kepler*np.cos(angle *np.pi/180)
        
        St     = (p.St)#.astype(str)
        r      = []
        t      = []
        Stlist = []

        r.append(p.r)
        t.append(p.t)
        Stlist.append(p.St)

        for idt, dtstep in enumerate(p.env.dt[:-1],1):
            while np.all(p.t < p.env.t[idt]):
                p.courant_velocity(same_dt = True)
                p.dt      = np.minimum(p.env.t[idt] - p.t, p.dt)
                p.dt[p.t >= p.env.t[idt]] = 0
                p.t = p.t + p.dt
                p.KDK_drag()
            p.env.next #NOTE det skal være uden for p.t comparision

            r.append(p.r)
            t.append(p.t)
            Stlist.append(p.St)
        
        r   = np.array(r); t = np.array(t); Stlist = np.array(Stlist)
        rr  = np.sqrt(np.sum(r**2,2))
        leg = ['%.2e'%(q) for q in St]
        # ------- r-t plots ----------
        ax = fig_rt.add_subplot(M, N ,idx+1)
        ax.set_title('angle = %i'%(angle))
            ## ---- Snap shot radius ----
        ax.plot(trc/runtime, rct/p.units.AU, label = 'rc' , color = 'grey', linewidth = 4)
            ## \------------------------/
        ax.plot(t[:,:,0]/runtime,rr/p.units.AU)
        ax.set_xlabel(r'time $[t_{collapse} = %.2e$ $yr]$'%(runtime/p.units.yr))                
        ax.set_ylabel(r'radius [AU]')
        ax.legend([r'$r_c(t)$']+leg,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_rt.tight_layout()
        # ------- r-r plots ------
        ax = fig_rr.add_subplot(M, N ,idx+1)
        ax.set_title('angle = %i'%(angle))
        ax.plot(r[:,:,0]/p.units.AU,r[:,:,1]/p.units.AU)
        ax.set_xlabel(r'x [AU]')                
        ax.set_ylabel(r'y [AU]')
        ax.set_xlim([-R/p.units.AU ,R/p.units.AU])
        ax.set_ylim([-R/p.units.AU ,R/p.units.AU])
        ax.legend(leg,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_rr.tight_layout()
        # ------- St-t plots ------#
        ax = fig_Stt.add_subplot(M,N, idx + 1)
        ax.set_title('angle = %i'%(angle))
        
        ax.plot(t[:,:,0]/runtime,Stlist[:,:,0])
        ax.set_xlabel(r'time $[t_{collapse} = %.2e$ $yr]$'%(runtime/p.units.yr))                
        ax.set_ylabel(r'St samples based on orbit radius')
        ax.legend(leg,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig_Stt.tight_layout()

        #------------------------
        print('angle %i done' %(angle))
    fig_rt.savefig('./Plots/collapse/%.2f-r-t.png'%(pr))
    fig_rr.savefig('./Plots/collapse/%.2f-r-r.png'%(pr))
    fig_Stt.savefig('./Plots/collapse/%.2f-St-t.png'%(pr))

if __name__ == '__main__':
    main()

# %% [markdown]
# # Legacy
# %% [markdown]
# def BE(units, M_BE = 1.0, T = 10):
#     """
#     Effectively a constuctor for the class. Most just so we dont have a long __init__, than can take up space
#     NOTE units is just sc.CGS look at the actuall scaling
#     TODO Calulate rc, rho_c and write then into the BE_sphere
#     TODO implement the use of the scale module.
#     TODO ! Have the lane-emden be save to a file, then by keyword or another function allow for regeneration of the values
#            Also save the crit values and allow the recalulation to find these again. This is not nesseary, but it is just nice to have, such that the file is lost
#            it can be regnerated. Recalulation every time is just a waste of time!!!
#     """
#     # check which M, R , T is none kinda like scaling does input check
#     mu = 2.42
#     Acrit,xcrit,ycrit,mcrit = np.loadtxt('./emden-solve/crit-value.txt',unpack = True) # Acrit, xcrit, ycrit, mcrit
#     xi,yi,_                 = np.loadtxt('./emden-solve/rho-adim.txt'  ,unpack = True)
#     M_BE  = M_BE * units.m_Sun
#     cs_sq = units.k_B*T/(mu*units.m_p)
#     rho_c = mcrit**2 * cs_sq**3/(4*np.pi*ycrit * M_BE**2 * units.G**3)
#     r_c   = np.sqrt(cs_sq/(4*np.pi * units.G * rho_c))
#     
#     sphere       = BE_sphere()
#     sphere.units = units
#     sphere.y     = scilate.interp1d(xi,yi,'linear')
#     sphere.xcrit = xcrit
#     sphere.xmin  = xi[0]; sphere.xmax  = xi[-1]
#     sphere.rc    = r_c
#     sphere.rho_c = rho_c
#     sphere.mu    = mu
#     return sphere
# %% [markdown]
# class BE_sphere(object):
#     # xmin, xmax are temp? maybe NOTE
#     __slots__ = 'units','rho_c','rc','xcrit','y', 'xmin','xmax','mu','rho','rho_d','T'
#                #sc.cgs , float, float, float, callable, float, float
#     def __init__(sphere,units, M_BE = 1.0, T = 10,drag=0.001):
#         Acrit,xcrit,ycrit,mcrit = np.loadtxt('./emden-solve/crit-value.txt',unpack = True) 
#         xi,yi,_                 = np.loadtxt('./emden-solve/rho-adim.txt'  ,unpack = True)
#         M_BE  = M_BE * units.m_Sun
#         cs_sq = units.k_B*T/(units.mu*units.m_p) # base units
#         rho_c = mcrit**2 * cs_sq**3/(4*np.pi*ycrit * M_BE**2 * units.G**3)
#         r_c   = np.sqrt(cs_sq/(4*np.pi * units.G * rho_c))
# 
#         sphere.T     = T
#         sphere.rho_d = 1.6# g/cm^3 # omregne TODO
#         sphere.units = units
#         sphere.rc    = r_c
#         sphere.rho_c = rho_c
#         sphere.mu    = units.mu
# 
#         sphere.y     = scilate.interp1d(xi,yi,'linear')
#         sphere.rho   = scilate.interp1d(xi*r_c,yi*rho_c,'linear')
#         sphere.xcrit = xcrit
#         sphere.xmin  = xi[0]; sphere.xmax  = xi[-1]
# 
#     def drag(self, r, s,rho_d = 1.6,):
#         """
#         rho_d = 1.6 g/cm^3
#         """
#         rm    = magnitude(r)
#         rho_g = self.rho(rm)
#         vth   = np.sqrt(8/np.pi) * self.units.k_B*self.T/(self.mu*self.units.m_p)
#         return rho_g*vth/(rho_d*s)
#         #vth =  self.units.k_B*self.T/(self.units.mu*self.units.m_p)
#  
# 
#     def radius_crit(self, N = 1000):
#         """
#         TODO: when scaling is properbly implemented, remove AU
#         """
#         x = np.linspace(self.xmin, self.xcrit, N)
#         return x*self.rc/self.units.AU
# 
#     def mass_total(self, N= 1000):
#         """
#         asssumes spherical symmetric in solar masses
#         TODO: when scaling is properbly implemented, remove m_Sun
#         """
#         x, dx = np.linspace(self.xmin, self.xcrit, N, retstep = True)
#         return 4*np.pi*np.cumsum(x**2*self.y(x))*dx* self.rc**3*self.rho_c/self.units.m_Sun       
# 
# 
# %% [markdown]
# class BE_sphere(object):
#     # xmin, xmax are temp? maybe NOTE
#     __slots__ = 'units','rho_c','rc','xcrit','y', 'xmin','xmax','mu','rho','T','M'
#                #sc.cgs , float, float, float, callable, float, float
#     def __init__(sphere,units, M_BE = 1.0, T = 10,drag=0.001):
#         Acrit,xcrit,ycrit,mcrit = np.loadtxt('./emden-solve/crit-value.txt',unpack = True) 
#         xi,yi,_                 = np.loadtxt('./emden-solve/rho-adim.txt'  ,unpack = True)
#         dx = xi[1]-xi[0]
# 
#         M_BE  = M_BE * units.m_Sun
#         cs_sq = units.k_B*T/(units.mu*units.m_p) # base units
#         rho_c = mcrit**2 * cs_sq**3/(4*np.pi*ycrit * M_BE**2 * units.G**3)
#         r_c   = np.sqrt(cs_sq/(4*np.pi * units.G * rho_c))
# 
#         
#         sphere.T     = T
#         sphere.units = units
#         sphere.rc    = r_c
#         sphere.rho_c = rho_c
#         sphere.mu    = units.mu
# 
#         m = 4*np.pi*np.cumsum(xi*xi * yi)*dx*rho_c * r_c**3
#         sphere.M     = scilate.interp1d(xi*r_c,m) ##TODO find a way to limit do xcrit?
# 
# 
#         sphere.y     = scilate.interp1d(xi,yi,'linear')
#         sphere.rho   = scilate.interp1d(xi*r_c,yi*rho_c,'linear')
#         sphere.xcrit = xcrit
#         sphere.xmin  = xi[0]; sphere.xmax  = xi[-1]
# 
#     def drag(self, r, s,rho_d = 1.6,):
#         """
#         rho_d = 1.6 g/cm^3
#         """
#         rm    = magnitude(r,True)
#         rho_g = self.rho(rm)
#         vth   = np.sqrt(8/np.pi * self.units.k_B*self.T/(self.mu*self.units.m_p))
#         return rho_g*vth/(rho_d*s)
#         #vth =  self.units.k_B*self.T/(self.units.mu*self.units.m_p)
# 
# 
# %% [markdown]
# class Particles(object):
#     #__slots__ = ''
#     # TODO maybe initalize everything as non ambigous vectors
#     def __init__(self, n, steps, env, smin, smax,rho_d = 1.6):
#         self.n  = n
#         self.t  = np.zeros((n,steps))
#         self.dt = np.zeros(self.n)
#         self.env = env ## enviorment ie sphere.
#         self.units = env.units
#         self.r = np.zeros((n,2))
#         self.v = self.r.copy()
#         self.s = np.linspace(smin,smax,n)
#         self.steps = steps
#         self.rho_d = rho_d # g/cm ^3
#         
#             
#     def courant_g(self,same_dt=True):
#         """ Determine time step(s) for at least self.nstep steps per orbit.
#             If same_dt is True, take the same time step for all particles,
#             If False, take nstep steps per orbit, for all radii
#         """
#         rm = magnitude(self.r)
#         vm = magnitude(self.v)
#         gm  = magnitude(self.gravity())
#         self.dt = 0.004 * vm/gm
#         #if same_dt:
#         #    self.dt[:] = self.dt.min()
#         #    self.t = self.t + self.dt[0]
# 
#     def gravity(self):
#         rm = magnitude(self.r)[:,None]
#         return  - self.r * self.units.G*self.units.m_Sun/(rm*rm*rm) # TODO add mass change
#         #return -self.r*self.G/(rm*rm*rm) # self .G is temp
# 
#     def kick(self,step=0.5):
#         """grav velocity step"""
#         force = self.gravity()#not really a force, more so an acceleration
#         self.v = self.v + step*self.dt[:,None]*force
#     
#     def drift(self,step=1.0):
#         """ Spacial step"""
#         self.r = self.r + step*self.dt*self.v
# 
#     def KDK(self,same_dt=True):
#         self.courant()
#         self.kick(0.5)
#         self.drift(1.0)
#         self.kick(0.5)
# 
#     def drag_debug(self, r, v):
#         """
#         returns the "drag velocity change" so \Delta v_{drag}
#         NOTE
#         Is here for debug reasons. The function call might not be worth it.
#         Although it should be minimal.
#         TODO
#         make this a function, that takes : r,p,v?
#         This does that we dont need to load this method as part of teh gas
#         """
#         Cdt    = self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
#         return Cdt[:,None]*(v - 0)/(1 + Cdt[:,None]) # vg = 0
# 
#     def KDK_drag(self,same_dt=True):
#         self.courant_g()
#         #self.dt = 3.154e+7
#         self.kick(0.5)
#         # NOTE + != * :(
#         r      = self.r  + 0.5 * self.dt * self.v # half here as it needs to sync with kick
#         Cdt    = self.dt * self.env.drag(r,self.s,self.rho_d) # NOTE changed self.r to r as it should be at time step 0.5
#         self.v = self.v  - Cdt*(self.v - 0)/(1 + Cdt) # vg = 0
#         self.drift(1.0)
#         self.kick(0.5)
#         Cdt    = self.dt * self.env.drag(self.r,self.s)
#         self.v = self.v  - Cdt*(self.v - 0)/(1 + Cdt) # vg = 0
# 
# %% [markdown]
# def main():
#     global p,v_kepler,r
#     sphere = BE_sphere(sc.CGS,1,10)
# 
#     p = Particles(1,3000,sphere,smin = 1e4, smax = 1e4)
# 
#     p.r[:,0] = 0.99*p.env.xcrit*p.env.rc
#     p.s[0]   = 1e4
#     angle    = 0
#     v_Kepler = np.sqrt(p.units.G*p.env.M(magnitude(p.r))/magnitude(p.r))
#     p.v[:,0] = -v_Kepler*np.sin(angle)
#     p.v[:,1] = +v_Kepler*np.cos(angle)
#     rs = magnitude(p.r)
# 
#     rv  = np.zeros((p.steps,2))
#     tv  = np.zeros((p.steps,1))
# 
#     for t in range(p.steps):
#         p.KDK_drag()
#         rv[t] = p.r
#         tv[t] = p.dt
# 
# 
#     with plt.style.context('ggplot'):
# 
#         fig = plt.figure(figsize = (19,9))
# 
#         ax = plt.subplot(131)
#         plt.title('s = %.2e cm'%(p.s[0]))
#         ax.axis('equal')
#         plt.xlim([-p.env.xcrit*p.env.rc/p.units.AU,p.env.xcrit*p.env.rc/p.units.AU])
#         plt.ylim([-p.env.xcrit*p.env.rc/p.units.AU,p.env.xcrit*p.env.rc/p.units.AU])
#         plt.plot(rv[:,0]/p.units.AU,rv[:,1]/p.units.AU)
# 
# 
#         ax = plt.subplot(132)
#         plt.plot(np.cumsum(tv), magnitude(rv))
#         plt.show()
# 
# if __name__ == '__main__':
#     main()
# %% [markdown]
# 

