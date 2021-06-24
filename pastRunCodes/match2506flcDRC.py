import numpy as np
import os

matchtol = 3

dir = 'catRawMags1305/catDir/'

cat = np.genfromtxt(dir+'flcDRCpos2506_f814w.dat')
drc = np.genfromtxt(dir+'drc_useful.dat')

RA, DEC, x, y, magr, id = 0, 1, 2, 3, 4, 5


RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4, xDRC, yDRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42


cat_DRCmags = np.zeros((len(cat),4))
cat_DRCmags[:,0] = np.arange(0,len(cat),1) # index, DRCx DRCy magDRC
cat = np.hstack((cat[:,[RA,DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, xDRC,yDRC]],cat_DRCmags))

RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, xDRC, yDRC, id_cat, xDRC_mat, yDRC_mat, magDRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15

master = drc[:,[x,y,magr]]

x,y,magr = 0,1,2

matchids = np.zeros((len(master),1)) #id_cat
#xDRC yDRC magDRC


nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xDRC]) <= matchtol) & (abs(master[row][y] - cat[:,yDRC]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_cat]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xDRC])**2 +  (master[row][y] - matchrows[dd][yDRC])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id_cat]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y magr id_cat'


np.savetxt(dir+'drcFLCidx_2506_f814w.dat',master,header=header)
