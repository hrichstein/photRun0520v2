import numpy as np

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
master = np.genfromtxt(dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT')
saveDir = './'

x, y, m606c, m814c, nstar, sat606, sat814, camera, m606, s606, q606, o606, f606, g606, rxs606, sky606, rmssky606, m814, s814, q814, o814, f814, g814, rxs814, sky814, rmssky814, ra, dec = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27

ra_arr = [43.8922300,43.9197460,43.9220320,43.8584350,43.8746860,43.8534360]
dec_arr = [-54.1360640,-54.1128260,-54.1165550,-54.1153960,-54.1152020,-54.1320310]

matchPix = np.zeros((len(master),6)) # array for the id,x,y values

matchtol = 1e-6

nF = True
row = int(0)

# print(master.shape)

for cc, coord in enumerate(ra_arr):
    print('CC=',cc)
    while (nF):
        matchrows = master[ (abs(master[row][ra]) - coord <= matchtol) & (abs(master[row][dec] - dec_arr[cc]) <= matchtol) ]

        if len(matchrows)==1:
            matchPix[row][0] =  matchrows[0][x]
            matchPix[row][1] =  matchrows[0][y]
            matchPix[row][2] =  matchrows[0][ra]
            matchPix[row][3] =  matchrows[0][dec]
            matchPix[row][4] =  coord
            matchPix[row][5] =  dec_arr[cc]
            row += 1
            print('Matched')

        elif len(matchrows) > 1:
            distDiff = np.zeros((len(matchrows),1))
            for mm in range(len(matchrows)):
                distDiff = np.sqrt( (master[row][ra] - matchrows[mm][ra])**2 + (master[row][dec] - matchrows[mm][dec])**2)
            small = np.argmin(distDiff)
            matchPix[row][0] = matchrows[small][x]
            matchPix[row][1] = matchrows[small][y]
            matchPix[row][2] =  matchrows[0][ra]
            matchPix[row][3] =  matchrows[0][dec]
            matchPix[row][4] =  coord
            matchPix[row][5] =  dec_arr[cc]
            row += 1

        else:
            matchPix = np.delete(matchPix, row, 0)
            master = np.delete(master, row, 0)
            row += 1
            # print('Deleted')

        if (row >= len(master)):
            nF = False

header = 'id x y ra dec ra_ref dec_ref'
form = '%d %1.6f %1.6f %1.7f %1.7f %1.7f %1.7f'

np.savetxt(saveDir+'pixMattiaWraDEC.dat',matchPix, fmt=form,header=header)
