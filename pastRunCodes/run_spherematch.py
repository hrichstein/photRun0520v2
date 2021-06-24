import numpy as np
from spherematch import *
import time

start=time.time()
targList = 'targnamesDirections.txt'
magDir = 'mags0811/magDir/'

names = np.genfromtxt(targList,dtype=str)
# names=['BOOTES-II-NORTH','BOOTES-II-SOUTH','CARINA','DRACO-II']

for nn in range(len(names)):
    infileI = names[nn]+'_F814W_magList.dat'
    infileV = names[nn]+'_F606W_magList.dat'

    xI = np.loadtxt(magDir+infileI, dtype='float')
    xV = np.loadtxt(magDir+infileV, dtype='float')

    raI = xI[:, 0]
    decI = xI[:, 1]
    raV = xV[:, 0]
    decV = xV[:, 1]

    idx1, idx2, ds = spherematch(raI, decI, raV, decV, tol=5e-5, nnearest=1)

    outfile = 'sp_match0811_3/'+names[nn]+'_F814W_match.txt'
    np.savetxt(outfile, idx1, fmt='%4i')
    outfile = 'sp_match0811_3/'+names[nn]+'_F606W_match.txt'
    np.savetxt(outfile, idx2, fmt='%4i')

print('It took {0:0.1f} seconds'.format(time.time() - start))
# print('----done !---')
