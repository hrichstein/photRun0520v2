import numpy as np

dir = 'catRawMags1305/catDir/'

master= np.genfromtxt(dir+'flcPSF_idx_1706_f814wAPER.dat')

cat = np.genfromtxt(dir+'flcPSFpos0906.dat')

# RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19, 20, 21

RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19

x,y,magr,id_cat = 0,1,2,3

newCols = np.zeros((len(master),14))

idCol = master[:,id_cat]
idx = np.asarray(idCol,int)

reg = cat[idx]

newCols[:,0] = reg[:,RA]
newCols[:,1] = reg[:,DEC]
newCols[:,2] = reg[:,flags]
newCols[:,3] = reg[:,c_star]
newCols[:,4] = reg[:,mag1]
newCols[:,5] = reg[:,mag2]
newCols[:,6] = reg[:,mag3]
newCols[:,7] = reg[:,mag4]
newCols[:,8] = reg[:,xt1]
newCols[:,9] = reg[:,yt1]
newCols[:,10] = reg[:,magZPT]
newCols[:,11] = reg[:,magZPTerr]
newCols[:,12] = reg[:,xPSF_trans]
newCols[:,13] = reg[:,yPSF_trans]


outArr = np.hstack((newCols,master[:,[x,y,magr,id_cat]]))
header= 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 magZPT magZPTerr xAPER_trans yAPER_trans xAPER_mas yAPER_mas magAPER id_cat'

np.savetxt(dir+'matchedFLCaper1706_814.dat',outArr,header=header)
