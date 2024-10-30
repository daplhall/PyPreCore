rp = np.loadtxt('../run067/output-00324.txt', usecols = (0,11,8,6,9))
cs = rp[:,4]
M = rp[:,2]
r = rp[:,0]
vkep = np.sqrt(G*M/r)
dm = M[1:] - M[:-1]
rho = rp[:,3]
dpdr = (rp[1:,1] - rp[:-1,1])/(rp[1:,0]-rp[:-1,0])


plt.figure(figsize = (14,6))

print(dpdr.max())
print(np.where(dpdr > 0), dpdr)

plt.plot(-dpdr)
plt.yscale('log')
plt.xlim([1700,1750])


dpdr = np.array([dpdr[0],dpdr])
np.array([dpdr[0],dpdr]).shape
dp = np.array([dpdr[0]].append(dpdr.tolist()))

dp = np.zeros(len(dpdr)+1)
dp[0] = dpdr[0] ; dp[1:] = dpdr

plt.figure(figsize = (14,6))
plt.plot(rho)
plt.yscale('log')

plt.figure(figsize = (14,6))
plt.plot(vkep)

plt.figure(figsize = (14,6))
plt.plot( 0.5*r*dp/(vkep*rho) ) ## Drift vel from 2010 bachelor article, eq 20, missing "-". Ed not used as it was a parameter for their model

plt.figure(figsize = (14,6))
plt.plot(cs**2/vkep*dp*rp[:,0]/rp[:,1]) # drift velocity

dpdr[0]
