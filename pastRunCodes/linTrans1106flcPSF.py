import numpy as np
from linear6d import *

dir = 'catRawMags1305/catDir/'


all_file = np.genfromtxt(dir+'flcAll1106_zptMags.dat')

RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1,yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18

#flc positions
x_arr = np.array([1942.41,2733.6358,3210.4213,275.484,1854.0268,3119.8681,3321.6462,1012.2739,2960.6373,1586.7523])

y_arr = np.array([ 336.0805,2591.7028,721.6712,2016.7489,2635.92,3469.8079,1810.6803,290.3904,3657.488,1715.3336])

match_arr = np.zeros((len(x_arr),2))

match_arr[:,0] = x_arr
match_arr[:,1] = y_arr

#psf positions
x_psf = np.array([2434.9197,4718.1953,2887.3555,4013.6179,4714.2852,5610.8672,3974.3972,2339.0757,5788.7881,3783.894])

y_psf = np.array([5232.4946,4569.6943,3994.6946,6979.1934,5445.4702,4233.6597,3943.2954,6153.9668,4402.1094,5660.8569])

master_arr = np.zeros((len(x_psf),2))

master_arr[:,0] = x_psf
master_arr[:,1] = y_psf

match = match_arr
master = master_arr
all = all_file

weights = np.zeros((len(match)))
weights.fill(1.0)

new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xt1],all[:,yt1])

# outname = 'flc2PSF_round3.dat'
outname = 'zptFlcPSF_1106.dat'

np.savetxt(dir+outname, new_all, fmt="%1.6f")

#################

cat = np.genfromtxt(dir+'flcAll1106_zptMags.dat')

transCat = np.genfromtxt(dir + 'zptFlcPSF_1106.dat')

newCol = np.zeros((len(cat),2))

newCol[:,0] = transCat[:,0]
newCol[:,1] = transCat[:,1]

cat = np.hstack((cat, newCol))

header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 mean stdev cut_flag idx_cut num_abv_std magZPT magZPTerr xPSF_trans yPSF_trans'

np.savetxt(dir+'flcPSFpos1106.dat', cat, header=header

# Iter      23    CHI-SQUARE =  2.52849533  DOF =  14
#    P0 = 1995.424216
#    P1 = 0.05447893527
#    P2 = 0.9930301159
#    P3 = 7143.23027
#    P4 = -0.9929655745
#    P5 = 0.05431881115
