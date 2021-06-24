import numpy as np
from scipy import stats

dir = 'catRawMags1305/catDir/'

# f814w = np.genfromtxt('flcDRC0506_magCut.dat',names=True)
# cat = np.genfromtxt('flcDRC0506_magCut.dat')
#
# all = np.genfromtxt('flcAll_magCut_0906.dat',names=True)
# all_cat = np.genfromtxt('flcAll_magCut_0906.dat')

f814w = np.genfromtxt('flcDRC1106_magCut.dat',names=True)
cat = np.genfromtxt('flcDRC1106_magCut.dat')

all = np.genfromtxt('flcAll_magCut_1106.dat',names=True)
all_cat = np.genfromtxt('flcAll_magCut_1106.dat')

RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16
# RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4, xDRC, yDRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42

# Calculate difference and corrections
flc_diff = stats.sigmaclip(f814w['magDRC']-f814w['mean'],2.5,2.5)
err_add = np.std(flc_diff[0] / np.sqrt(len(flc_diff[0])))
mag_corr = np.nanmean(flc_diff[0])

# Apply to all
new_mag = all['mean'] + mag_corr
new_err = np.sqrt( all['stdev']**2 + err_add**2 )

newCol = np.zeros((len(all),2))
newCol[:,0] = new_mag
newCol[:,1] = new_err

out = np.hstack((all_cat,newCol))

header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 mean stdev cut_flag idx_cut num_abv_std magZPT magZPTerr'

np.savetxt(dir+'flcAll1106_zptMags.dat',out,header=header)
