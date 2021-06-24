import numpy as np
import matplotlib.pyplot as plt
import os
from linear6d import *

def linFiltTransDRC(targname,dir='./'):

    x = os.listdir(dir)
    for ii in x:
        if ii.endswith('_F606WphotU.dat'):
            f606w_file = ii
        elif ii.endswith('_F814WphotU.dat'):
            f814w_file = ii

    # Going from F814W to F606W
    file = np.genfromtxt(dir+'drcFiltRef_'+targname+'_2.dat',names=True)
    fileCat = np.genfromtxt(dir+'drcFiltRef_'+targname+'_2.dat')
    colNs = np.array(file.dtype.names)

    all = np.genfromtxt(dir+f814w_file,names=True)
    allCat = np.genfromtxt(dir+f814w_file)
    colAs = np.array(all.dtype.names)

    # Putting the F814W positions into the match array to be used as references transformed
    match_arr = np.zeros((len(file),2))
    xi_2v = np.int(np.where(colNs=='x_f814w')[0])
    yi_2v = np.int(np.where(colNs=='y_f814w')[0])

    match_arr[:,0] = fileCat[:,xi_2v]
    match_arr[:,1] = fileCat[:,yi_2v]

    # Putting the F606W positions into the master array
    master_arr = np.zeros((len(file),2))
    xv = np.int(np.where(colNs=='x_f606w')[0])
    yv = np.int(np.where(colNs=='y_f606w')[0])

    master_arr[:,0] = fileCat[:,xv]
    master_arr[:,1] = fileCat[:,yv]

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    # Filling the array with values to be transformed
    all_arr = np.zeros((len(all),2))

    x_bt = np.int(np.where(colAs=='xcenter')[0])
    y_bt = np.int(np.where(colAs=='ycenter')[0])

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    outName = dir + targname + "_DRCfiltTrans2"

    new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

    makePlot(targname,match_arr[:,0],match_arr[:,1],\
    new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in F814W',label_2='New in F814W 2 F606W',label_3='Original in F606W',outname=outName+'_matchCheck')

    outArr = np.hstack((allCat,new_all))

    s0=' '
    header = s0.join(colAs)
    header += ' x_f606wTrans y_f606wTrans'

    # np.savetxt(outName+'.dat',outArr,header=header)
    np.savetxt(outName+'_pU.dat',outArr,header=header)

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
    plt.savefig(outname+'_pU.png',dpi=600,bbox_inches='tight')
    plt.close()


    return None
