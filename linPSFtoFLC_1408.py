import numpy as np
import matplotlib.pyplot as plt
from linear6d import *

def linPSFTrans(targname,dir='photUtils0820/'):

    # Going from PSF to FLC_pu

    file = np.genfromtxt(dir+'flcPSFref_'+targname+'_1408.dat',names=True)
    fileCat = np.genfromtxt(dir+'flcPSFref_'+targname+'_1408.dat')
    colNs = np.array(file.dtype.names)

    dir2 = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
    all = np.genfromtxt(dir2+'HOROLOGIUM_CF.1.TOSEND.CAT',names=True)
    allCat = np.genfromtxt(dir2+'HOROLOGIUM_CF.1.TOSEND.CAT')
    colAs = np.array(all.dtype.names)

    # Putting the PSF positions into the match array to be used as references transformed
    match_arr = np.zeros((len(file),2))
    xa_2p = np.int(np.where(colNs=='xPSF')[0])
    ya_2p = np.int(np.where(colNs=='yPSF')[0])

    match_arr[:,0] = fileCat[:,xa_2p]
    match_arr[:,1] = fileCat[:,ya_2p]

    # Putting the PU positions into the master array
    master_arr = np.zeros((len(file),2))
    xp = np.int(np.where(colNs=='x_f606w')[0])
    yp = np.int(np.where(colNs=='y_f606w')[0])

    master_arr[:,0] = fileCat[:,xp]
    master_arr[:,1] = fileCat[:,yp]

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    # Filling the array with values to be transformed
    all_arr = np.zeros((len(all),2))

    x_bt = np.int(np.where(colAs=='x')[0])
    y_bt = np.int(np.where(colAs=='y')[0])

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    outName = dir + targname + "_PSF2flcTrans"

    new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])


    outArr = np.hstack((allCat,new_all))

    s0=' '
    header = s0.join(colAs)
    header += ' x_Trans y_Trans'

    # np.savetxt(outName+'.dat',outArr,header=header)
    np.savetxt(outName+'_1408.dat',outArr,header=header)

    makePlot(targname,match_arr[:,0],match_arr[:,1],\
    new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in PSF',label_2='New in PSF 2 FLC',label_3='Original in FLC',outname=outName+'_matchCheck')

    return None


def makePlot(targname,x1,y1,x2,y2,x3,y3,label_1,\
    label_2,label_3,outname):

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(x3,y3,label=label_3,s=70)
    # ax.scatter(x1,y1,label=label_1,s=50)
    ax.scatter(x2,y2,label=label_2,s=20)

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outname+'_pU.png',dpi=600,bbox_inches='tight')
    plt.close()


    return None
