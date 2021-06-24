import numpy as np
import os

matchtol = 1.5

dir = 'catRawMags1305/catDir/'

psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT')

aper = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.2.TOSEND.CAT')

x, y, m606c, m814c, nstar, sat606, sat814, camera, m606, s606, q606,o606, f606, g606,rxs606,sky606,rmssky606, m814,s814 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 #psf

id_aper = np.zeros((len(aper),1))
id_aper[:,0] = np.arange(0,len(aper),1)

aper = np.hstack((aper[:,[x,y,nstar]],id_aper))

master = psf[:,[x, y, m606c, m814c,s606,s814,nstar]]

x, y, nstarAPER, idAPER = 0,1,2,3

cat = aper

matchids = np.zeros((len(master),2)) #id_cat

nF = True
row = 0

while (nF):
    matchrows = cat[(abs(master[row][x] - cat[:,x]) <= matchtol) & (abs(master[row][y] - cat[:,y]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][nstarAPER]
      matchids[row][1] = matchrows[0][idAPER]
      # match_info[row][0] = master[row][x]
      # match_info[row][1] = master[row][y]
      # match_info[row][2] = master[row][magr]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][x])**2 +  (master[row][y] - matchrows[dd][y])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][nstarAPER]
         matchids[row][1] = matchrows[small][idAPER]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)
        # match_info = np.delete(match_info,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x y m606c m184c s606 s814 nstarPSF nstarAPER idAPER'

#
np.savetxt(dir+'psfAPER_idx_1706.dat',master,header=header)
