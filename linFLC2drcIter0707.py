import numpy as np
from linear6d import *
import matplotlib.pyplot as plt

def linFLC2drc_i(targname,filt,dir='./',iter=1):

    if iter==1:

        # file = np.genfromtxt(dir+'flcDRCref_'+filt+".dat",names=True)
        # fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+".dat")

        # all = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans.dat",names=True)
        # allCat = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans.dat")
        file = np.genfromtxt(dir+'flcDRCref_'+filt+"_mDc.dat",names=True)
        fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+"_mDc.dat")

        all = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans_mDc.dat",names=True)
        allCat = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans_mDc.dat")

        colAs = np.array(all.dtype.names)

        # outName = dir+targname+"_"+filt+"_drcTrans"+str(iter)
        # plotName = outName+'_matchCheck_'+str(iter)

        s0=' '
        header = s0.join(colAs)
        header += ' xDRC_{0:d} yDRC_{0:d}'.format(iter)

        x_bt = np.int(np.where(colAs=='xDRC')[0])
        y_bt = np.int(np.where(colAs=='yDRC')[0])

    elif iter> 1:

        # file = np.genfromtxt(dir+'flcDRCref_'+filt+"_"+str(iter-1)+'.dat',names=True)
        # fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+"_"+str(iter-1)+'.dat')
        #
        # all = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter-1)+".dat",names=True)
        # allCat = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter-1)+".dat")

        file = np.genfromtxt(dir+'flcDRCref_'+filt+"_"+str(iter-1)+'_mDc.dat',names=True)
        fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+"_"+str(iter-1)+'_mDc.dat')

        all = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter-1)+"_mDc.dat",names=True)
        allCat = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter-1)+"_mDc.dat")

        colAs = np.array(all.dtype.names)

        # outName = dir+targname+"_"+filt+"_drcTrans"+str(iter)
        # plotName = outName+'_matchCheck_'+str(iter)

        s0=' '
        header = s0.join(colAs)
        header += ' xDRC_{0:d} yDRC_{0:d}'.format(iter)

        itx = 'xDRC_'+ str(iter-1)
        ity = 'yDRC_'+ str(iter-1)

        x_bt = np.int(np.where(colAs==itx)[0])
        y_bt = np.int(np.where(colAs==ity)[0])

    outName = dir+targname+"_"+filt+"_drcTrans"+str(iter)
    plotName = outName+'_matchCheck_'+str(iter)

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

    all_arr = np.zeros((len(all),2))

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    try:
        new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

        # This would be a good place to have it make a plot

        makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
        new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Past FLC Trans Pos.',label_2='New in FLC2DRC',label_3='Original in DRC',outname=plotName)

        outArr = np.hstack((allCat,new_all))
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

#    # if iter==None:
    #     file = np.genfromtxt(dir+'flcDRCref_'+filt+'.dat',names=True)
    #     fileCat = np.genfromtxt(dir+'flcDRCref_'+filt+'.dat')
    #
    #     all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    #     allCat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')
    #
    #     colAs = np.array(all.dtype.names)
    #     x_bt = np.int(np.where(colAs=='xt1')[0])
    #     y_bt = np.int(np.where(colAs=='yt1')[0])
    #
    #     outName = dir+targname+"_"+filt+"_drcTrans"
    #     plotName = outName+'_matchCheck.png'
    #
    #     s0=' '
    #     header = s0.join(colAs)
    #     header += ' xDRC yDRC'
