import numpy as np

dir = 'catRawMags1305/catDir/'

# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'flcPSF_idx_2506_all.dat')

flc = np.genfromtxt(dir+'flcPSFpos2506_all.dat') # 606
psf = np.genfromtxt(dir+'matchedPSFaper1706_tc.dat') #814

x,y,id_flc,id_psf = 0,1,2,3

idCol_flc = master[:,id_flc]
idx_flc = np.asarray(idCol_flc,int)
reg_flc = flc[idx_flc]

##################################

# newCols814 = np.zeros((len(master),11))
idCol_psf = master[:,id_psf]
idx_psf = np.asarray(idCol_psf,int)
reg_psf = psf[idx_psf]

outArr = np.hstack((reg_flc,reg_psf))
header= 'c_star_f606w mag1_f606w mag2_f606w mag3_f606w mag4_f606w xr1_f606w yr1_f606w xt1_f606w yt1_f606w xDRC_trans_f606w yDRC_trans_f606w xDRC_mat_f606w yDRC_mat_f606w magDRC_f606w id_cat_f606w mean_f606w stdev_f606w cut_flag_f606w idx_cut_f606w num_abv_std_f606w magZPT_f606w magZPTerr_f606w c_star_f814w mag1_f814w mag2_f814w mag3_f814w mag4_f814w xr1_f814w yr1_f814w xt1_f814w yt1_f814w xDRC_trans_f814w yDRC_trans_f814w xDRC_mat_f814w yDRC_mat_f814w magDRC_f814w id_cat_f814w mean_f814w stdev_f814w cut_flag_f814w idx_cut_f814w num_abv_std_f814w magZPT_f814w magZPTerr_f814w xPSF_trans yPSF_trans xPSF yPSF m606cPSF m814cPSF s606PSF s814PSF nstarPSF nstarAPER idAPER xAPER yAPER m606cAPER m814cAPER s606APER  s814APER '

np.savetxt(dir+'matchedFLCpsf2506_all.dat',outArr,header=header)
