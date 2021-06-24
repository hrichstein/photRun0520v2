import numpy as np
from linear6d import *

dir = 'catRawMags1305/catDir/'


all_file = np.genfromtxt(dir+'flcDRC0506_zptMags.dat')

RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, xDRC_trans, yDRC_trans, xDRC_mat, yDRC_mat, magDRC, id_cat, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
#flc positions

x_arr = np.array([2733.5338,1735.3688,1586.7124,1870.5519,1578.4385,1843.006,1942.4388])

y_arr = np.array([2591.8732,1915.1304,1715.1404	,1625.8403,736.27513	,625.11832,336.7416])

match_arr = np.zeros((len(x_arr),2))

match_arr[:,0] = x_arr
match_arr[:,1] = y_arr

#psf positions

x_psf = np.array([4718.1953,3991.6089,3783.894,3711.5825,2812.5161,2716.2554,2434.9197 ])
y_psf = np.array([4569.6943	,5523.9956	,5660.8569	,5374.0845	,5615.8149	,5347.2842	,5232.4946])

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
outname = 'zptFlcPSF.dat'

np.savetxt(dir+outname, new_all, fmt="%1.6f")

#################

cat = np.genfromtxt(dir+'flcDRC0506_zptMags.dat')

transCat = np.genfromtxt(dir + 'zptFlcPSF.dat')

newCol = np.zeros((len(cat),2))

newCol[:,0] = transCat[:,0]
newCol[:,1] = transCat[:,1]

cat = np.hstack((cat, newCol))

header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 xDRC_trans yDRC_trans xDRC_mat yDRC_mat magDRC id_cat mean stdev cut_flag idx_cut num_abv_std magZPT xPSF_trans yPSF_trans'

np.savetxt(dir+'flcPSFpos0506.dat', cat, header=header)
