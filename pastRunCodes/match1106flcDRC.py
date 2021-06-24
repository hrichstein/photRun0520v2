import numpy as np
import os

matchtol = 3

dir = 'catRawMags1305/catDir/'

# dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
# psf = np.genfromtxt(dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT')

# cat = np.genfromtxt(dir+'flcDRCpos1106_f606w.dat')
cat = np.genfromtxt(dir+'flcDRCpos1106_f606wr2.dat')

drc = np.genfromtxt(dir+'drc_useful_f606w.dat')

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
# match_info = np.zeros((len(master),3))
# cat_pos = np.zeros((len(master),2))

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xDRC]) <= matchtol) & (abs(master[row][y] - cat[:,yDRC]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_cat]
      # match_info[row][0] = master[row][x]
      # match_info[row][1] = master[row][y]
      # match_info[row][2] = master[row][magr]
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
        # match_info = np.delete(match_info,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))

# print(len(cat))
# matchids= np.asarray(matchids,dtype=int)
# for ii, idx in enumerate(matchids):
#     cat[idx,xDRC_mat] = match_info[ii][x]
#     cat[idx,yDRC_mat] = match_info[ii][y]
#     cat[idx,magDRC] = match_info[ii][magr]
#
#     cat[idx,xDRC_mat] = match_info[ii][x]
#     cat[idx,yDRC_mat] = match_info[ii][y]
#     cat[idx,magDRC] = match_info[ii][magr]
#
#
# cat = cat[cat[:,magr]>=1e-3]
# print(len(cat))

master= np.hstack((master,matchids))
header= 'x y magr id_cat'

# cat[:,xDRC_mat] = match_info[:,x]
# cat[:,yDRC_mat] = match_info[:,y]
# cat[:,magDRC] = match_info[:,magr]
#
#
# header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 xDRC yDRC id_cat xDRC_mat yDRC_mat magDRC'
#
np.savetxt(dir+'drcFLCidx_1106_f606wr2.dat',master,header=header)
