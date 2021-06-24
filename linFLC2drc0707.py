import numpy as np
from linear6d import *
import matplotlib.pyplot as plt

def linFLC2drc(targname,filt,dir='./'):

    # file = np.genfromtxt(dir+'flcDRCref_'+filt+'.dat',names=True)
    # fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+'.dat')

    file = np.genfromtxt(dir+'flcDRCref_'+filt+'_mDc.dat',names=True)
    fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+'_mDc.dat')

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

    # all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    # allCat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')

    all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'_mDc.dat',names=True)
    allCat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'_mDc.dat')

    colAs = np.array(all.dtype.names)
    x_bt = np.int(np.where(colAs=='xt1')[0])
    y_bt = np.int(np.where(colAs=='yt1')[0])

    all_arr = np.zeros((len(all),2))

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    outName = dir+targname+"_"+filt+"_drcTrans"

    try:
        new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

        # This would be a good place to have it make a plot

        # makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
        # new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in FLC',label_2='New in FLC2DRC',label_3='Original in DRC',outname=outName+'_matchCheck.png')

        makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
        new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in FLC',label_2='New in FLC2DRC',label_3='Original in DRC',outname=outName+'_matchCheck_mDc')

        outArr = np.hstack((allCat,new_all))
        header = 'RA DEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xr2 yr2 xr3 yr3 xr4 yr4 xc1 yc1 xc2 yc2 xc3 yc3 xc4 yc4 xt1 yt1 xt2 yt2 xt3 yt3 xt4 yt4 mean stdev cut_flag idx_cut num_abv_std xDRC yDRC'
        #
        # np.savetxt(outName+'.dat',outArr,header=header)
        np.savetxt(outName+'_mDc.dat',outArr,header=header)

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

    plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
    plt.close()


    return None

#
