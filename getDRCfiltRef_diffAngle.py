import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits


def getRefDRCFilt(targname,dir='./',matchtol=5):

    infoFile = np.genfromtxt('../Hannah_Data/drcTargInfo_new.dat',names=True)

    x = os.listdir(dir)
    for ii in x:
        if ii.endswith('_F606WphotU.dat'):
            f606w_file = ii
        elif ii.endswith('_F814WphotU.dat'):
            f814w_file = ii

    f606wN = np.genfromtxt(dir+f606w_file,names=True)
    f606w = np.genfromtxt(dir+f606w_file)

    f814wN = np.genfromtxt(dir+f814w_file,names=True)
    f814w = np.genfromtxt(dir+f814w_file)

    col606 = np.array(f606wN.dtype.names)
    col814 = np.array(f814wN.dtype.names)

    # Getting indices of columns

    # Will be the same for both filters
    x = np.int(np.where(col606=='xcenter')[0])
    y = np.int(np.where(col606=='ycenter')[0])
    mag = np.int(np.where(col606=='magr')[0])

    idCol = np.int(np.where(col606=='id')[0])

    # Sorting to get the 50 brightest stars
    v50 = np.argsort(f606wN['magr'])[:50]
    fv_50 = f606w[v50]

    i50 = np.argsort(f814wN['magr'])[:50]
    fi_50 = f814w[i50]

    master_in = fv_50[:,[idCol,x,y,mag]]
    idV,xv,yv,magv = 0,1,2,3

    cat = fi_50
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,xv,yv,x,y,idCol)

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = fv_50[:,[idCol,x,y,mag]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 3

    master = np.hstack((master,matchids))

    idV, xv, yv, magv, idI = 0, 1, 2, 3, 4

    newCols = np.zeros((len(master),3))

    idxCol = master[:,idI]
    idxI = np.asarray(idxCol,int)
    regI = f814w[idxI]

    newCols[:,0] = regI[:,x]
    newCols[:,1] = regI[:,y]
    newCols[:,2] = regI[:,mag]

    outArr = np.hstack((master,newCols))

    idV, xv, yv, magv, idI, xi, yi, magi = 0, 1, 2, 3, 4, 5, 6, 7

    header = 'id_f606w x_f606w y_f606w magr_f606w id_f814w x_f814w y_f814w magr_f814w'
    form = '%d %1.5f %1.5f %1.4f %d %1.5f %1.5f %1.4f'

    outName = dir+'drcFiltRef_'+targname
    # np.savetxt(outName+'.dat',outArr,header=header,fmt=form)
    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xv],outArr[:,yv],label='F606W',s= 50)
    ax.scatter(outArr[:,xi],outArr[:,yi],label='F814W',s=20)

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.close()


    return None


def matchlistID(master,cat,matchtol,x1,y1,x2,y2,id_mat):

    matchids_in = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:

        matchrows = cat[(abs(master[row][x1] - cat[:,x2]) \
            <= matchtol) & (abs(master[row][y1] - cat[:,y2])<= matchtol)]

        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[0][id_mat]
            row += 1

        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt( (master[row][x1] - \
                matchrows[dd][x2])**2 +  (master[row][y1] \
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

                # print(len(udx),len(master))
                master = master[udx]
                matchids_in = matchids_in[udx]

                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(matchtol,len(master)))
                nF = False

            elif len(udx)==len(master):
                # print(len(udx),len(master))
                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(matchtol,len(master)))
                nF = False

    return master,matchids_in


#
