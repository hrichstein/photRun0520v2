import numpy as np

# dir = 'catRawMags1305/catDir/'

dir = 'sDRC_2606/'
# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'drcFullidx_2606.dat')


f606w = np.genfromtxt(dir+'drcDRCpos2606_intoF814W.dat')
f814w = np.genfromtxt(dir+'cat_HOR-I_F814W_cStarCut.dat')

# There should be a last column in each catalog giving the index placement

x,y,mag814_in,id814,id606 = 0,1,2,3,4
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
header= 'flags_f606w RA_f606w DEC_f606w xr_f606w yr_f606w flux_f606w c_star_f606w magr_f606w id_f606w xr_f814w_trans yr_f814w_trans flags_f814w RA_f814w DEC_f814w xr_f814w yr_f814w flux_f814w c_star_f814w magr_f814w id_f814w'

np.savetxt(dir+'matchedDRCfullCat_2606.dat',outArr,header=header,fmt='%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d')
