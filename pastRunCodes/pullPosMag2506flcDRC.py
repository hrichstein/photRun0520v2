import numpy as np

dir = 'catRawMags1305/catDir/'

# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'drcDRC_idx_2506_mas.dat')


f606w = np.genfromtxt(dir+'flcAll2506_zptMags_f606w.dat') # 606
f814w = np.genfromtxt(dir+'flcAll2506_zptMags_f814w.dat') #814

# There should be a last column in each catalog giving the index placement

c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, xDRC_trans, yDRC_trans, xDRC_mat, yDRC_mat, magDRC, id_cat, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22

x,y,mag606_in,id606,id814 = 0,1,2,3,4
# Master header
# x y magZ id606 id814

# newCols606 = np.zeros((len(master),11))

idCol606 = master[:,id606]
idx606 = np.asarray(idCol606,int)
reg606 = f606w[idx606]

##################################

# newCols814 = np.zeros((len(master),11))
idCol814 = master[:,id814]
idx814 = np.asarray(idCol814,int)
reg814 = f814w[idx814]

outArr = np.hstack((reg606,reg814))
header= 'c_star_f606w mag1_f606w mag2_f606w mag3_f606w mag4_f606w xr1_f606w yr1_f606w xt1_f606w yt1_f606w xDRC_trans_f606w yDRC_trans_f606w xDRC_mat_f606w yDRC_mat_f606w magDRC_f606w id_cat_f606w mean_f606w stdev_f606w cut_flag_f606w idx_cut_f606w num_abv_std_f606w magZPT_f606w magZPTerr_f606w c_star_f814w mag1_f814w mag2_f814w mag3_f814w mag4_f814w xr1_f814w yr1_f814w xt1_f814w yt1_f814w xDRC_trans_f814w yDRC_trans_f814w xDRC_mat_f814w yDRC_mat_f814w magDRC_f814w id_cat_f814w mean_f814w stdev_f814w cut_flag_f814w idx_cut_f814w num_abv_std_f814w magZPT_f814w magZPTerr_f814w'

np.savetxt(dir+'matchedFLCdrc2506_comb.dat',outArr,header=header)
