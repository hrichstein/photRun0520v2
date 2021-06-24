import numpy as np

dir3  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
drc = np.genfromtxt(dir3 + 'HOROLOGIUM-I_sfErr.dat',names=True)

idx= np.argsort(drc['magRaw_i'])[:1000]

drc_50 = drc[idx]


new_cols = np.zeros((len(drc_50),2))
new_cols[:,0] = drc_50['x_i']
new_cols[:,1] = drc_50['y_i']

np.savetxt('drc1000.reg',new_cols,fmt='%1.5f')
