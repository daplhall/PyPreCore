import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
import scipy.optimize as sciopt
import scipy.interpolate as scilate
import seaborn as sns
x = np.array([1,2,3,4,5,6,7,8, 2,9,0,76,5])
xx = np.linspace(0,10, len(x))
f = scilate.interp1d(xx,x,copy = False)
x[0] = -50

f(xx[0])