import numpy as np
import os

matchtol = 2

dir = 'catRawMags1305/catDir/'

aper_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
aper_file = np.genfromtxt(aper_dir + 'HOROLOGIUM_CF.2.TOSEND.CAT')

cat = np.genfromtxt(dir+'drc2APER_f814w_t.dat')

cat_id = np.zeros((len(cat),1))
cat_id[:,0] = np.arange(0,len(cat),1)

cat = np.hstack((cat,cat_id))

x, y, m606c, m814c = 0, 1, 2, 3 #aper

RA_v, DEC_v, x_v, y_v, fAper_v, fErr_v, magAper_v, magErr_v, magRaw_v, magRed_v, magAbs_v, elong_v, ellip_v, class_Star_v, RA_i, DEC_i, x_i, y_i, fAper_i, fErr_i, magAper_i, magErr_i, magRaw_i, magRed_i, magAbs_i, elong_i, ellip_i, class_Star_i, corrF_errV, corrF_errI, corrM_errV, corrM_errI,xt_i,yt_i,id = 0, 1, 2 ,3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34


master = aper_file[:,[x,y,m814c]]

x,y,magr = 0,1,2

matchids = np.zeros((len(master),1)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,xt_i]) <= matchtol) & (abs(master[row][y] - cat[:,yt_i]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][id]

      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xt_i])**2 +  (master[row][y] - matchrows[dd][yt_i])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][id]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)
        # match_info = np.delete(match_info,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y magr_814 id_cat'

# cat[:,xDRC_mat] = match_info[:,x]
# cat[:,yDRC_mat] = match_info[:,y]
# cat[:,magDRC] = match_info[:,magr]
#
#
# header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xt1 yt1 xDRC yDRC id_cat xDRC_mat yDRC_mat magDRC'
#
np.savetxt(dir+'drcAPER_idx_2406_f814w.dat',master,header=header)
