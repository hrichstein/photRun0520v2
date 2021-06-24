import numpy as np

dir = 'catRawMags1305/catDir/'

# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'flcFLC_idx_1206_OC.dat')

f606w = np.genfromtxt(dir+'matchedFLCpsf1106.dat') # 606
f814w = np.genfromtxt(dir+'matchedFLCpsf0906.dat') #814

# There should be a last column in each catalog giving the index placement


# F814W
RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1_814, yt1_814, magZPT814,magZPTerr814, xPSF_trans814, yPSF_trans814, xPSF_mas814, yPSF_mas814, magPSF814, id_cat = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17

# RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 magZPT magZPTerr xPSF_trans yPSF_trans xPSF_mas yPSF_mas magPSF id_cat

#F606W
RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1_606, yt1_606, magZPT606, magZPTerr606, xPSF_trans606, yPSF_trans606, xPSF_mas606, yPSF_mas606, magPSF606, id_new606 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19

x,y,mag606,id606,id814 = 0,1,2,3,4
# Master header
# x y magZ id606 id814

newCols606 = np.zeros((len(master),17))

idCol606 = master[:,id606]
idx606 = np.asarray(idCol606,int)
reg606 = f606w[idx606]

newCols606[:,0] = reg606[:,RA]
newCols606[:,1] = reg606[:,DEC]
newCols606[:,2] = reg606[:,flags]
newCols606[:,3] = reg606[:,c_star]
newCols606[:,4] = reg606[:,mag1]
newCols606[:,5] = reg606[:,mag2]
newCols606[:,6] = reg606[:,mag3]
newCols606[:,7] = reg606[:,mag4]
newCols606[:,8] = reg606[:,xr1]
newCols606[:,9] = reg606[:,yr1]
newCols606[:,10] = reg606[:,xt1_606]
newCols606[:,11] = reg606[:,yt1_606]
newCols606[:,12] = reg606[:,magZPT606]
newCols606[:,13] = reg606[:,magZPTerr606]
newCols606[:,14] = reg606[:,xPSF_trans606]
newCols606[:,15] = reg606[:,yPSF_trans606]
newCols606[:,16] = reg606[:,magPSF606]
##################################

newCols814 = np.zeros((len(master),15))
idCol814 = master[:,id814]
idx814 = np.asarray(idCol814,int)
reg814 = f814w[idx814]

newCols814[:,0] = reg814[:,RA]
newCols814[:,1] = reg814[:,DEC]
newCols814[:,2] = reg814[:,flags]
newCols814[:,3] = reg814[:,c_star]
newCols814[:,4] = reg814[:,mag1]
newCols814[:,5] = reg814[:,mag2]
newCols814[:,6] = reg814[:,mag3]
newCols814[:,7] = reg814[:,mag4]
newCols814[:,8] = reg814[:,xt1_814]
newCols814[:,9] = reg814[:,yt1_814]
newCols814[:,10] = reg814[:,magZPT814]
newCols814[:,11] = reg814[:,magZPTerr814]
newCols814[:,12] = reg814[:,xPSF_trans814]
newCols814[:,13] = reg814[:,yPSF_trans814]
newCols814[:,14] = reg814[:,magPSF814]

outArr = np.hstack((newCols606,newCols814))
header= 'RA_f606w DEC_f606w flags_f606w c_star_f606w mag1_f606w mag2_f606w mag3_f606w mag4_f606w xr1_f606w yr1_f606w xt1_f606w yt1_f606w magZPT_f606w magZPTerr_f606w xPSF_trans_f606w yPSF_trans_f606w magPSF_f606w RA_f814w DEC_f814w flags_f814w c_star_f814w mag1_f814w mag2_f814w mag3_f814w mag4_f814w xt1_f814w yt1_f814w magZPT_f814w magZPTerr_f814w xPSF_trans_f814w yPSF_trans_f814w magPSF_f814w'

# np.savetxt(dir+'matchedPSFpsf1106.dat',outArr,header=header)
np.savetxt(dir+'matchedFLCflc1206.dat',outArr,header=header)
