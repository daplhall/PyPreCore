# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import Simulation as sim
import scaling as sc
sb.set(context = 'notebook',font = 'serif',font_scale = 1.3)
sb.set_palette('Set2', 8)


# %%
ginterface  = np.geomspace(1.42921456e-06, 1.42921456e-04, 20)


# %%

path = './SimData/Fullrun_10/'


# %%
sphere = sim.DataSphere   ('./run_extracted/run067/', units = sc.CGS)
sphere._vg[:,:-1] = (sphere._rc[:,1:]-sphere._rc[:,:-1])/sphere.dt[None,:-1]
sphere.update()


# %%
rowrange = { 'skiprows' : 0,    # Start
             'max_rows' : 30 } # \# Lines to read after start
p        = sim.DataParticles(path, sphere, **rowrange)

# %%
D = sim.dust_to_gas_DataParticles(p,p.s[0].min(),29)




