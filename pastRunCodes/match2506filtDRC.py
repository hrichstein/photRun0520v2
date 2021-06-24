import numpy as np
import os

matchtol = 3

dir = 'catRawMags1305/catDir/'


f606w = np.genfromtxt(dir+'flcAll2506_zptMags_f606w.dat') # 606
f814w = np.genfromtxt(dir+'flcAll2506_zptMags_f814w.dat') #814

cat_id = np.zeros((len(f814w),1))
cat_id[:,0] = np.arange(0,len(f814w),1)
#
f814w = np.hstack((f814w,cat_id))

cat_id2 = np.zeros((len(f606w),1))
cat_id2[:,0] = np.arange(0,len(f606w),1)
f606w = np.hstack((f606w,cat_id2))

# F814W
c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, xDRC_trans, yDRC_trans, xDRC_mat, yDRC_mat, magDRC, id_cat, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22


#F606W
c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, xDRC_trans, yDRC_trans, xDRC_mat, yDRC_mat, magDRC, id_cat, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22


master = f606w[:,[xDRC_mat, yDRC_mat, magZPT, id_new]]
x,y,magr = 0,1,2

cat = f814w

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xDRC_mat]) <= matchtol) & (abs(master[row][y] - cat[:,yDRC_mat]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_new]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xDRC_mat])**2 +  (master[row][y] - matchrows[dd][yDRC_mat])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id_new]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y magZ_f606w id606 id814'

#
np.savetxt(dir+'drcDRC_idx_2506_mas.dat',master,header=header)
