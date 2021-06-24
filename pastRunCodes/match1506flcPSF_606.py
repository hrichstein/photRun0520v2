import numpy as np
import os

matchtol = 2

dir = 'catRawMags1305/catDir/'

psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf_file = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT')

cat = np.genfromtxt(dir+'flcPSFpos1106.dat')

cat_id = np.zeros((len(cat),1))
cat_id[:,0] = np.arange(0,len(cat),1)

cat = np.hstack((cat,cat_id))

x, y, m606c, m814c, nstar, sat606, sat814, camera, m606, s606, q606,o606, f606, g606,rxs606,sky606,rmssky606, m814,s814 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 #psf

# flc
RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19, 20, 21


master = psf_file[:,[x,y,m606c,s606,nstar]]

x,y,magr,magerr,n_id = 0,1,2,3,4

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xPSF_trans]) <= matchtol) & (abs(master[row][y] - cat[:,yPSF_trans]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_new]
      # match_info[row][0] = master[row][x]
      # match_info[row][1] = master[row][y]
      # match_info[row][2] = master[row][magr]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xPSF_trans])**2 +  (master[row][y] - matchrows[dd][yPSF_trans])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id_new]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)
        # match_info = np.delete(match_info,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y magr magerr nstar id_cat'

# cat[:,xDRC_mat] = match_info[:,x]
# cat[:,yDRC_mat] = match_info[:,y]
# cat[:,magDRC] = match_info[:,magr]
#
#
# header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 xDRC yDRC id_cat xDRC_mat yDRC_mat magDRC'
#
np.savetxt(dir+'flcPSF_idx_1506_606.dat',master,header=header)
