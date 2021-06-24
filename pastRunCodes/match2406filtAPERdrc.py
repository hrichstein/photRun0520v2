import numpy as np
import os

matchtol = 2

dir = 'catRawMags1305/catDir/'

# psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
# psf_file = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT')

# f606w = np.genfromtxt(dir+'matchedFLCpsf1506_606.dat') # 606
# f814w = np.genfromtxt(dir+'matchedFLCpsf1506_814.dat') #814

f606w = np.genfromtxt(dir+'matchedDRCaper2406_606.dat') # 606
f814w = np.genfromtxt(dir+'matchedDRCaper2406_814.dat') #814

cat_id = np.zeros((len(f814w),1))
cat_id[:,0] = np.arange(0,len(f814w),1)
#
f814w = np.hstack((f814w,cat_id))

cat_id2 = np.zeros((len(f606w),1))
cat_id2[:,0] = np.arange(0,len(f606w),1)
f606w = np.hstack((f606w,cat_id2))

mag606, magErr606, xt_v, yt_v, xAPER_606, yAPER_606, magAPER_606, idAPER_606,id_new606 = 0, 1, 2, 3, 4, 5, 6, 7,8

mag814, magErr814, xt_v, yt_v, xAPER_814, yAPER_814, magAPER_814, idAPER_814,id_new814 = 0, 1, 2, 3, 4, 5, 6, 7,8


master = f606w[:,[xAPER_606, yAPER_606, mag606, id_new606]]
x,y,magr = 0,1,2

cat = f814w

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xAPER_814]) <= matchtol) & (abs(master[row][y] - cat[:,yAPER_814]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id_new814]
      # match_info[row][0] = master[row][x]
      # match_info[row][1] = master[row][y]
      # match_info[row][2] = master[row][magr]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xAPER_814])**2 +  (master[row][y] - matchrows[dd][yAPER_814])**2)
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
header= 'x y mag606 id606 id814'


np.savetxt(dir+'drcAPER_idx_2406_mas.dat',master,header=header)
