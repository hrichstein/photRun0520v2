from astropy.io import fits
import numpy as np
import os

dir = 'catRawMags1305/catDir/'

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf = np.genfromtxt(dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT')

def matchXYs(master,match,matchtol=1,outName='welp'):

    idx_mat = np.zeros((len(match), 1))
    idx_mat[:,0] = np.arange(0,len(match),1)

    idx_mas = np.zeros((len(master), 1))
    idx_mas[:,0] = np.arange(0,len(master),1)

    mas_ids = np.zeros((len(master), 1))

    matcher = np.hstack((match[:,[0,1]],idx_mat))
    # matcher = np.hstack((match[:,[0,1]],idx_mat))

    # matchCat = matcher.copy()

    nF = True
    row = 0

    while (nF):
        matchrows = match[ (abs(master[row][0] - match[:,0]) <= matchtol)  & (abs(master[row][1] - match[:,1]) <= matchtol) ]

        # matchrows = match[ (abs(y_mas[row] - y_mat[row]) <= matchtol) ]

        if (len(matchrows) == 1):
            mas_ids[row] = matcher[row][2]
            row += 1
            # matchids[row] = matcher[row][2]
            # row = row + 1
            # print('match')

        # elif (len(matchrows) > 1):
        #     distDiff = np.zeros((len(matchrows),1))
        #     for dd in range(len(matchrows)):
        #         distDiff[dd] = np.sqrt( (master[row][0] - match[:,0])**2 +  (master[row][1] - match[:,1])**2)
        #     small = np.argmin(distDiff)
        #     matchids[row] = idx_mas[small]
        #     row += 1

        else:
            print('delete')
            print(row)
            # print(len(matcher))
            master = np.delete(master, row, 0)
            mas_ids = np.delete(mas_ids,row,0)
            matcher = np.delete(matcher,row,0)
            match = np.delete(match, row, 0)
            # matchids = np.delete(matchids,row,0)
            # matchCat = np.delete(matchCat,row,0)
            idx_mat = np.delete(idx_mat,row,0)

        if (row >= len(master)):
            nF = False
            print('Master Length:',len(master))

    outArr = np.hstack((master,mas_ids,matcher))

    np.savetxt(dir+outName+'.dat',outArr,fmt='%1.5f %1.5f %d', header='x y matchID')

    return None

master = psf[:,[0,1]]
match = np.genfromtxt(dir+'flc2PSF_round1.dat')

matchXYs(master,match,matchtol=1,outName='flc2PSF_match_r1_0406')
