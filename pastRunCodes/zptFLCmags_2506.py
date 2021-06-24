import numpy as np
from scipy import stats

dir = 'catRawMags1305/catDir/'


f606w = np.genfromtxt('flcDRC2506_magCut_f606w.dat',names=True)
cat = np.genfromtxt('flcDRC2506_magCut_f606w.dat')

c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, xDRC_trans, yDRC_trans, xDRC_mat, yDRC_mat, magDRC, id_cat, mean, stdev, cut_flag, idx_cut, num_abv_std = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19


# Calculate difference and corrections
flc_diff = stats.sigmaclip(f606w['magDRC']-f606w['mean'],2.5,2.5)
err_add = np.std(flc_diff[0] / np.sqrt(len(flc_diff[0])))
mag_corr = np.nanmean(flc_diff[0])

# Apply to all
new_mag = f606w['mean'] + mag_corr
new_err = np.sqrt( f606w['stdev']**2 + err_add**2 )

newCol = np.zeros((len(f606w),2))
newCol[:,0] = new_mag
newCol[:,1] = new_err

out = np.hstack((cat,newCol))

header = 'c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 xDRC_trans yDRC_trans xDRC_mat yDRC_mat magDRC id_cat mean stdev cut_flag idx_cut num_abv_std magZPT magZPTerr'

np.savetxt(dir+'flcAll2506_zptMags_f606w.dat',out,header=header)
