import numpy as np
from astropy.stats import sigma_clip

dir= 'catRawMags1305/catDir/'

file = np.genfromtxt(dir+'matched_w_MagsPos2705r3.dat')


RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40


mag_arr = np.zeros((len(file),4))
mag_arr[:,0] = file[:,mag1]
mag_arr[:,1] = file[:,mag2]
mag_arr[:,2] = file[:,mag3]
mag_arr[:,3] = file[:,mag4]

mag_median = np.nanmedian(mag_arr,axis=1)

newCol = np.zeros((len(file),5))

newCol[:,0] = file[:,xt1]
newCol[:,1] = file[:,yt1]
newCol[:,2] = mag_median
newCol[:,3] = file[:,RA]
newCol[:,4] = file[:,DEC]

np.savetxt(dir+'pixPosHorImattia.dat',newCol,header='x y mag ra dec',fmt='%1.6f %1.6f %1.5f %1.7f %1.7f')
