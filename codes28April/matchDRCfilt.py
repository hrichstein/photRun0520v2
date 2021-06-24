"""
Matching sources between DRC filters using 6D-transformed positions
In the previous code, I transformed F814W sources into F606W space.
I will take the original F606W catalog and the new F814W catalog.
"""

import numpy as np
import os


def matchFiltDRC(targname,dir='./',matchtol=2,suffix='_pu.dat'):

    x = os.listdir(dir)
    for ii in x:
        if ii.endswith('_F606W' + suffix):
            f606w_file = ii
        # elif ii.endswith('_F814W' + suffix):
        #     f814w_file = ii

    f606wN = np.genfromtxt(dir + f606w_file,names=True)
    f606w = np.genfromtxt(dir + f606w_file)

    f814wN = np.genfromtxt(dir + "drcFiltTrans_" + targname + ".dat",
                           names=True)
    f814w = np.genfromtxt(dir + "drcFiltTrans_" + targname + ".dat")

    col606 = np.array(f606wN.dtype.names)
    col814 = np.array(f814wN.dtype.names)

    # Getting columns of x,y of transformed F814W to F606W
    xI = np.int(np.where(col814=='x_f606wTrans')[0])
    yI = np.int(np.where(col814=='y_f606wTrans')[0])
    idI = np.int(np.where(col814=='id')[0])

    xV = np.int(np.where(col606=='xcenter')[0])
    yV = np.int(np.where(col606=='ycenter')[0])
    idV = np.int(np.where(col606=='id')[0])

    master_in = f606w[:,[idV,xV,yV]]
    idV_m,x,y = 0,1,2

    cat = f814w
    nF_out = True

    len606 = len(f606w)
    len814 = len(f814w)

    minLen = np.min([len606,len814])

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xI,yI,idI)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,
                                                             len(master)))
            matchtol += 1
            if matchtol < 7:
                master_in = f606w[:,[idV,xV,yV]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    print(len(master)/minLen)  # the percentage of sources matched
    master = np.hstack((master,matchids))

    idV_mas, xV_mas, yV_mas, idI_mas = 0, 1, 2, 3

    idCol606 = master[:,idV_mas]
    idx606 = np.asarray(idCol606,int)
    reg606 = f606w[idx606]

    idCol814 = master[:,idI_mas]
    idx814 = np.asarray(idCol814,int)
    reg814 = f814w[idx814]

    outArr = np.hstack((reg606,reg814))

    s6 = '_f606w '
    header606 = s6.join(col606)

    s8 = '_f814w '
    header814 = s8.join(col814)

    header = header606 + '_f606w ' + header814

    np.savetxt(dir + targname + '_matchedDRCfilt.dat',outArr,header=header)

    return None


def matchlistID(master,cat,matchtol,x1,y1,x2,y2,id_mat):

    matchids_in = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:

        matchrows = cat[(abs(master[row][x1] - cat[:,x2])
                        <= matchtol) & (abs(master[row][y1] - cat[:,y2])
                        <= matchtol)]

        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[0][id_mat]
            row += 1

        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt((master[row][x1]
                                       - matchrows[dd][x2])**2
                                       + (master[row][y1]
                                       - matchrows[dd][y2])**2)
            small = np.argmin(distDiff)
            matchids_in[row][0] = matchrows[small][id_mat]
            row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]
                nF = False

            elif len(udx)==len(master):
                nF = False

    return master,matchids_in

#
