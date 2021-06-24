import numpy as np
import os

matchtol = 2

dir = 'catRawMags1305/catDir/'

# psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
# psf_file = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT')

f606w = np.genfromtxt(dir+'matchedFLCpsf1106.dat') # 606
f814w = np.genfromtxt(dir+'matchedFLCpsf0906.dat') #814

cat_id = np.zeros((len(f814w),1))
cat_id[:,0] = np.arange(0,len(f814w),1)
#
f814w = np.hstack((f814w,cat_id))

cat_id2 = np.zeros((len(f606w),1))
cat_id2[:,0] = np.arange(0,len(f606w),1)
f606w = np.hstack((f606w,cat_id2))

# cat = np.hstack((cat,cat_id))
# F814W
RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, magZPT814,magZPTerr, xPSF_trans814, yPSF_trans814, xPSF_mas814, yPSF_mas814, magPSF, id_cat, id_new814 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18

# RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 magZPT magZPTerr xPSF_trans yPSF_trans xPSF_mas yPSF_mas magPSF id_cat

#F606W
RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, magZPT606, magZPTerr606, xPSF_trans606, yPSF_trans606, xPSF_mas606, yPSF_mas606, magPSF,id_cat, id_new606 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19, 20

# RA DEC flags c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 magZPT magZPTerr xPSF_trans yPSF_trans xPSF_mas yPSF_mas magPSF id_cat
# RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19, 20, 21


master = f606w[:,[xPSF_mas606, yPSF_mas606, magZPT606, id_new606]]
x,y,magr = 0,1,2

cat = f814w

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xPSF_mas814]) <= matchtol) & (abs(master[row][y] - cat[:,yPSF_mas814]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_new814]
      # match_info[row][0] = master[row][x]
      # match_info[row][1] = master[row][y]
      # match_info[row][2] = master[row][magr]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xPSF_mas814])**2 +  (master[row][y] - matchrows[dd][yPSF_mas814])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id_new814]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)
        # match_info = np.delete(match_info,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y magZ id606 id814'

# cat[:,xDRC_mat] = match_info[:,x]
# cat[:,yDRC_mat] = match_info[:,y]
# cat[:,magDRC] = match_info[:,magr]
#
#
# header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 xDRC yDRC id_cat xDRC_mat yDRC_mat magDRC'
#
np.savetxt(dir+'psfPSF_idx_1106_mas.dat',master,header=header)
