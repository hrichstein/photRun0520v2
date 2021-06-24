import numpy as np
from linear6d import *
import matplotlib.pyplot as plt

def linFLC2drc(targname,filt,dir='./'):

    file = np.genfromtxt(dir+'flcDRCref_'+filt+'_pu.dat',names=True)
    fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+'_pu.dat')

    colNs = np.array(file.dtype.names)

    match_arr = np.zeros((len(file),2))
    xt = np.int(np.where(colNs=='xF')[0])
    yt = np.int(np.where(colNs=='yF')[0])

    match_arr[:,0] = fileCat[:,xt]
    match_arr[:,1] = fileCat[:,yt]

    master_arr = np.zeros((len(file),2))
    xt1 = np.int(np.where(colNs=='xD')[0])
    yt1 = np.int(np.where(colNs=='yD')[0])

    master_arr[:,0] = fileCat[:,xt1]
    master_arr[:,1] = fileCat[:,yt1]

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    allCat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')


    colAs = np.array(all.dtype.names)

    s0 = ' '
    header = s0.join(colAs)

    x_bt = np.int(np.where(colAs=='xt1')[0])
    y_bt = np.int(np.where(colAs=='yt1')[0])

    all_arr = np.zeros((len(all),2))

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    outName = dir+targname+"_"+filt+"_drcTrans"

    try:
        new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

        makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
        new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in FLC',label_2='New in FLC2DRC',label_3='Original in DRC',outname=outName+'_matchCheck')

        outArr = np.hstack((allCat,new_all))
        header += ' xDRC yDRC'

        np.savetxt(outName+'_pu.dat',outArr,header=header)

    except RuntimeWarning:
        print('Not good enough.',targname,filt,dd)


    return None


def makePlot(targname,filt,x1,y1,x2,y2,x3,y3,label_1,\
    label_2,label_3,outname):

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(x3,y3,label=label_3,s=70)
    ax.scatter(x1,y1,label=label_1,s=50)
    ax.scatter(x2,y2,label=label_2,s=20)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    plt.savefig(outname+'_pu.png',dpi=600,bbox_inches='tight')
    plt.close()


    return None

#
