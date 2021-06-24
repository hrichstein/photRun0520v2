import numpy as np
import os

matchtol = 3

dir = 'catRawMags1305/catDir/'


flc = np.genfromtxt(dir+'flcPSFpos2506_all.dat')
psf = np.genfromtxt(dir+'matchedPSFaper1706_tc.dat')

cat_id = np.zeros((len(flc),1))
cat_id[:,0] = np.arange(0,len(flc),1)
#
flc = np.hstack((flc,cat_id))

cat_id2 = np.zeros((len(psf),1))
cat_id2[:,0] = np.arange(0,len(psf),1)
psf = np.hstack((psf,cat_id2))

# flc
c_star_f606w, mag1_f606w, mag2_f606w, mag3_f606w, mag4_f606w, xr1_f606w, yr1_f606w, xt1_f606w, yt1_f606w, xDRC_trans_f606w, yDRC_trans_f606w, xDRC_mat_f606w, yDRC_mat_f606w, magDRC_f606w, id_cat_f606w, mean_f606w, stdev_f606w, cut_flag_f606w, idx_cut_f606w, num_abv_std_f606w, magZPT_f606w, magZPTerr_f606w, c_star_f814w, mag1_f814w, mag2_f814w, mag3_f814w, mag4_f814w, xr1_f814w, yr1_f814w, xt1_f814w, yt1_f814w, xDRC_trans_f814w, yDRC_trans_f814w, xDRC_mat_f814w, yDRC_mat_f814w, magDRC_f814w, id_cat_f814w, mean_f814w, stdev_f814w, cut_flag_f814w, idx_cut_f814w, num_abv_std_f814w, magZPT_f814w, magZPTerr_f814w, xPSF, yPSF, id_flc = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46

# psf
xPSF_mat,        yPSF_mat,        m606cPSF,   m814cPSF,   s606PSF,       s814PSF,     nstarPSF,   nstarAPER,   idAPER,   xAPER,       yAPER,       m606cAPER,   m814cAPER,   s606APER,    s814APER, id_psf = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15


master = flc[:,[xPSF, yPSF, id_flc]] # shorter one
x,y = 0,1

cat = psf

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xPSF_mat]) <= matchtol) & (abs(master[row][y] - cat[:,yPSF_mat]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_psf]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xPSF_mat])**2 +  (master[row][y] - matchrows[dd][yPSF_mat])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id_psf]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'xPSF_trans yPSF_trans id_flc id_psf'

#
np.savetxt(dir+'flcPSF_idx_2506_all.dat',master,header=header)
