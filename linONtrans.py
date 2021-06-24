import numpy as np
import matplotlib.pyplot as plt
from linear6d import *

def linONtrans(targname,filt='F606W',dir='./'):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    # Going from Old to New
    file = np.genfromtxt(dir+'oldNewref_'+filt+'.dat',names=True)
    fileCat = np.genfromtxt(dir+'oldNewref_'+filt+'.dat')
    colNs = np.array(file.dtype.names)

    allN = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I/HOROLOGIUM-I_allMatchedZPTed.dat',names=True)
    allCat = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I/HOROLOGIUM-I_allMatchedZPTed.dat')
    colOs = np.array(allN.dtype.names)

    # Putting the F814W positions into the match array to be used as references transformed
    match_arr = np.zeros((len(file),2))
    xo_2n = np.int(np.where(colNs=='xO')[0])
    yo_2n = np.int(np.where(colNs=='yO')[0])

    match_arr[:,0] = fileCat[:,xo_2n]
    match_arr[:,1] = fileCat[:,yo_2n]

    # Putting the F606W positions into the master array
    master_arr = np.zeros((len(file),2))
    xv = np.int(np.where(colNs=='xN')[0])
    yv = np.int(np.where(colNs=='yN')[0])

    master_arr[:,0] = fileCat[:,xv]
    master_arr[:,1] = fileCat[:,yv]

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    # Filling the array with values to be transformed
    all_arr = np.zeros((len(allN),2))

    oxstr = 'xt1'+fils
    oystr = 'yt1'+fils
    x_bt = np.int(np.where(colOs==oxstr)[0])
    y_bt = np.int(np.where(colOs==oystr)[0])

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    outName = dir + targname + "_ONtrans_"+filt

    new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

    makePlot(targname,match_arr[:,0],match_arr[:,1],\
    new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in Old',label_2='New in Old 2 New',label_3='Original in New',outname=outName+'_matchCheck')

    outArr = np.hstack((allCat,new_all))

    s0=' '
    header = s0.join(colOs)
    header += ' xt1_ON yt1_ON'

    # np.savetxt(outName+'.dat',outArr,header=header)
    np.savetxt(outName+'.dat',outArr,header=header)

    return None


def makePlot(targname,x1,y1,x2,y2,x3,y3,label_1,\
    label_2,label_3,outname):

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(x3,y3,label=label_3,s=70)
    ax.scatter(x1,y1,label=label_1,s=50)
    ax.scatter(x2,y2,label=label_2,s=20)

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
    plt.close()


    return None
