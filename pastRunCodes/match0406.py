import numpy as np
import os

matchtol = 2.8

dir = 'catRawMags1305/catDir/'

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf = np.genfromtxt(dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT')

# cat = np.genfromtxt(dir+'flc2PSF_round1.dat')
# cat = np.genfromtxt(dir+'flc2PSF_round2.dat')
cat = np.genfromtxt(dir+'flc2PSF_round4.dat')
cat_ids = np.zeros((len(cat),1))
cat_ids[:,0] = np.arange(0,len(cat),1)
cat = np.hstack((cat,cat_ids))

master = psf[:,[0,1]]

cat_ids_m = np.zeros((len(master),1))
mas_ids = np.zeros((len(master),1))
mas_ids[:,0] = np.arange(0,len(master),1)

master = np.hstack((master,mas_ids,cat_ids_m))

cat_pos = np.zeros((len(master),2))

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][0] - cat[:,0]) <= matchtol) & (abs(master[row][1] - cat[:,1]) <= matchtol)]

    if (len(matchrows) == 1):
      master[row][3] = matchrows[0][2] #putting in id of star in cat
      cat_pos[row][0] = matchrows[0][0]
      cat_pos[row][1] = matchrows[0][1]
      row = row + 1

    else:
        master = np.delete(master,row,0)
        cat_pos = np.delete(cat_pos,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))

master = np.hstack((master,cat_pos))

header = 'x_psf y_psf id_psf id_cat x_cat y_cat'
form = '%1.5f %1.5f %d %d %1.5f %1.5f'

# np.savetxt(dir+'match_0406.dat',master,header=header,fmt=form)

# np.savetxt(dir+'match_0406r2.dat',master,header=header,fmt=form)
# 669

# np.savetxt(dir+'match_0406r4.dat',master,header=header,fmt=form)
# 669
