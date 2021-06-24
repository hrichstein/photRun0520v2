""" Getting references to transform from FLC to DRC and transforming """

import numpy as np

import matplotlib.pyplot as plt

from linear6d import test_linear
from matchlistID import matchlistID


def transFLCdrc(targname,flcFile,drcFile,dir='./',matchtol=1):

    flcN = np.genfromtxt(dir+flcFile,names=True)
    flc = np.genfromtxt(dir+flcFile)

    drcN = np.genfromtxt(dir+drcFile,names=True)
    drc = np.genfromtxt(dir+drcFile)

    colFs = np.array(flcN.dtype.names)
    colDs = np.array(drcN.dtype.names)

    xF = np.int(np.where(colFs=='xDRC1_f606w')[0])
    yF = np.int(np.where(colFs=='yDRC1_f606w')[0])

    xD = np.int(np.where(colDs=='xcenter_f606w')[0])
    yD = np.int(np.where(colDs=='ycenter_f606w')[0])

    idColF = len(colFs)
    newCol = np.zeros((len(flc),1),dtype=int)
    newCol[:,0] = np.arange(0,len(flc),1)

    flc_id = np.hstack((flc,newCol))

    idColD = len(colDs)
    newCol = np.zeros((len(drc),1),dtype=int)
    newCol[:,0] = np.arange(0,len(drc),1)
    drc_id = np.hstack((drc,newCol))

    f50 = np.argsort(flcN['mean_f606w'])[:50]
    flc_50 = flc_id[f50]

    d50 = np.argsort(drcN['magr_f606w'])[:50]
    drc_50 = drc_id[d50]

    master_in = drc_50[:,[idColD,xD,yD]]

    idD, xd, yd = 0, 1, 2

    cat = flc_50

    nF_out = True

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,xd,yd,xF,yF,
                                       idColF)

        if len(master)>=int(6):  # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = drc_50[:,[idColD,xD,yD]]  # resetting, just in case
            matchtol += 1

    master = np.hstack((master,matchids))
    idD, xd, yd, idF = 0, 1, 2, 3

    newCols = np.zeros((len(master),2))
    idxCol = master[:,idF]
    idxF = np.asarray(idxCol,int)
    regF = flc[idxF]

    newCols[:,0] = regF[:,xF]
    newCols[:,1] = regF[:,yF]

    tempArr = np.hstack((master,newCols))

    idD, xd, yd, idF, xf, yf = 0, 1, 2, 3, 4, 5

    # Linear Transform Part

    match_arr = np.zeros((len(tempArr),2))
    match_arr[:,0] = tempArr[:,xf]  # the x in FLC that is going to DRC
    match_arr[:,1] = tempArr[:,yf]  # the y in FLC that is going to DRC

    master_arr = np.zeros((len(tempArr),2))
    master_arr[:,0] = tempArr[:,xd]  # the x ref DRC
    master_arr[:,1] = tempArr[:,yd]  # the y ref DRC

    weights = np.ones(len(master_arr))

    all_arr = np.zeros((len(flc),2))
    all_arr[:,0] = flc_id[:,xF]
    all_arr[:,1] = flc_id[:,yF]

    print('Transforming ',targname)
    new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1],
                                     master_arr[:,0],master_arr[:,1],weights,
                                     weights, all_arr[:,0],all_arr[:,1])

    outArr = np.hstack((flc,new_all))
    s0 = ' '
    header = s0.join(colFs)
    header += ' x_DRCtrans y_DRCtrans'

    outName = dir + 'flcDRCtrans_' + targname

    np.savetxt(outName + '.dat',outArr,header=header)

    makePlot(targname,match_arr[:,0],match_arr[:,1],master_arr[:,0],
             master_arr[:,1],new_match[:,0],new_match[:,1],
             label_1='Original in FLC',label_2='Original in DRC',
             label_3='New in FLC 2 DRC',outname=outName+'_matchCheck')

    return None


def makePlot(targname,x1,y1,x2,y2,x3,y3,label_1,
             label_2,label_3,outname=None):

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(x1,y1,label=label_1,s=60)
    ax.scatter(x2,y2,label=label_2,s=25)
    ax.scatter(x3,y3,label=label_3,s=10)

    ax.legend()
    ax.set_title(targname)

    plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
    plt.close()

    return None

#
