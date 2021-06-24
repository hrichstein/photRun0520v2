import numpy as np

matchtol = 3

dir = 'sDRC_2606/'
dir2 = 'catRawMags1305/catDir/'
# dir2 = 'catRawMags1305/catDir/'

drc = np.genfromtxt(dir+'aperDRCpos2606_fullTransAPER.dat')
aper = np.genfromtxt(dir2+'matchedPSFaper1706_tc.dat')

drc_id = np.zeros((len(drc),1))
drc_id[:,0] = np.arange(0,len(drc),1) # index, DRCx DRCy magDRC
drc = np.hstack((drc,drc_id))

flags_f606w, RA_f606w, DEC_f606w, xr_f606w, yr_f606w, flux_f606w, c_star_f606w, magr_f606w, id_f606w, xr_f814w_trans, yr_f814w_trans, flags_f814w, RA_f814w, DEC_f814w, xr_f814w, yr_f814w, flux_f814w, c_star_f814w, magr_f814w, id_f814w, xAPER_trans, yAPER_trans, idNew_DRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22

aper_id = np.zeros((len(aper),1))
aper_id[:,0] = np.arange(0,len(aper),1) # index, DRCx DRCy magDRC
aper = np.hstack((aper,aper_id))

# APER heading

xPSF, yPSF, m606cPSF, m814cPSF, s606PSF, s814PSF, nstarPSF, nstarAPER, idAPER, xAPER, yAPER, m606cAPER, m814cAPER, s606APER, s814APER, idNew_APER =  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15

master = aper[:,[xAPER,yAPER,idNew_APER]]

x,y = 0,1

matchids = np.zeros((len(master),1)) #id_f606w
#xDRC yDRC magDRC

nF = True
row = 0

while (nF):
    matchrows = drc[(abs(master[row][x] - drc[:,xAPER_trans]) <= matchtol) & (abs(master[row][y] - drc[:,yAPER_trans]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][idNew_DRC]
      row = row + 1

    elif (len(matchrows) > 1):
         distDiff = np.zeros((len(matchrows),1))
         for dd in range(len(matchrows)):
             distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xAPER_trans])**2 +  (master[row][y] - matchrows[dd][yAPER_trans])**2)
         small = np.argmin(distDiff)
         matchids[row][0] = matchrows[small][idNew_DRC]
         row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x_APER y_APER id_APER id_DRC'


np.savetxt(dir+'aperDRCidx_2606_all.dat',master,header=header)
