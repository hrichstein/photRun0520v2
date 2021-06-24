import numpy as np

# dir = 'catRawMags1305/catDir/'

dir = 'sDRC_2606/'
dir2 = 'catRawMags1305/catDir/'
# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'aperDRCidx_2906_all_bkgd.dat')

drc = np.genfromtxt(dir+'aperDRCpos2906_fullTransAPER_bkgd.dat')
aper = np.genfromtxt(dir2+'matchedPSFaper2906_tc.dat')

x_APER, y_APER, id_aper, id_drc = 0, 1, 2, 3
# Master header
# x y magZ id606 id814

# newCols606 = np.zeros((len(master),11))

idCol_aper = master[:,id_aper]
idx_aper = np.asarray(idCol_aper,int)
reg_aper = aper[idx_aper]

##################################

# newCols814 = np.zeros((len(master),11))
idCol_drc = master[:,id_drc]
idx_drc = np.asarray(idCol_drc,int)
reg_drc = drc[idx_drc]

outArr = np.hstack((reg_aper,reg_drc))
header= 'xPSF yPSF m606cPSF m814cPSF s606PSF s814PSF nstarPSF nstarAPER idAPER xAPER yAPER m606cAPER m814cAPER s606APER s814APER sky606APER   sky814APER ra dec flags_f606w RA_f606w DEC_f606w xr_f606w yr_f606w flux_f606w bkgd_f606w c_star_f606w magr_f606w id_f606w xr_f814w_trans yr_f814w_trans flags_f814w RA_f814w DEC_f814w xr_f814w yr_f814w flux_f814w bkgd_f814w _star_f814w magr_f814w id_f814w xAPER_trans yAPER_trans'

np.savetxt(dir+'matchedDRCaper2906_bkgd.dat',outArr,header=header)
