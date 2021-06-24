import numpy as np

file = np.genfromtxt('flc6source.dat')

RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40


np.savetxt('flcXY.dat',file[:,[RA,DEC,xr1,yr1,xt1,yt1]],header='RA DEC xr1 yr1 xt1 yt1',fmt='%1.7f %1.7f %1.5f %1.5f %1.5f %1.5f')
