import numpy as np

dir = 'catRawMags1305/catDir/'

master= np.genfromtxt(dir+'drcFLCidx_2506_f814w.dat')

x,y,magr,id_cat = 0,1,2,3

cat = np.genfromtxt(dir+'flcDRCpos2506_f814w.dat')

RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4, xDRC, yDRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42

newCols = np.zeros((len(master),11))

idCol = master[:,id_cat]
idx = np.asarray(idCol,int)

reg = cat[idx]

newCols[:,0] = reg[:,c_star]
newCols[:,1] = reg[:,mag1]
newCols[:,2] = reg[:,mag2]
newCols[:,3] = reg[:,mag3]
newCols[:,4] = reg[:,mag4]
newCols[:,5] = reg[:,xr1]
newCols[:,6] = reg[:,yr1]
newCols[:,7] = reg[:,xt1]
newCols[:,8] = reg[:,yt1]
newCols[:,9] = reg[:,xDRC]
newCols[:,10] = reg[:,yDRC]

outArr = np.hstack((newCols,master[:,[x,y,magr,id_cat]]))
header= 'c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 xDRC_trans yDRC_trans xDRC_mat yDRC_mat magDRC id_cat'

np.savetxt(dir+'matchedFLCdrc2506_f814w.dat',outArr,header=header)
