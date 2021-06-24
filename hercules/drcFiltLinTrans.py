""" Code to do 6D linear transform """

import numpy as np
import matplotlib.pyplot as plt
import os

from linear6d import test_linear


def drcFiltLinTrans(targname,f8d,dir='./',suffix='_pu.dat',setnum=0):

    # Will need to replace filter names with appropriate
    # ones for your use.
    x = os.listdir(dir)
    for ii in x:
        # if ii.endswith('_F606W' + suffix):
        #     f606w_file = ii
        if ii.endswith(f8d+'_drc_F814W' + suffix):
            f814w_file = ii

    # Going from F814W to F606W
    refN = np.genfromtxt(dir + 'drcFiltRef_' + targname + str(setnum) + '.dat',
                         names=True)
    ref = np.genfromtxt(dir + 'drcFiltRef_' + targname + str(setnum) + '.dat')
    colRs = np.array(refN.dtype.names)

    allN = np.genfromtxt(dir + f814w_file,names=True)
    all = np.genfromtxt(dir + f814w_file)
    colAs = np.array(allN.dtype.names)

    # Putting the F814W positions into the match array
    # to be used as references.
    match_arr = np.zeros((len(ref),2))
    xi2v = np.int(np.where(colRs=='x_f814w')[0])
    yi2v = np.int(np.where(colRs=='y_f814w')[0])

    match_arr[:,0] = ref[:,xi2v]
    match_arr[:,1] = ref[:,yi2v]

    # Putting the F606W positions into the master array
    master_arr = np.zeros((len(ref),2))
    xv = np.int(np.where(colRs=='x_f606w')[0])
    yv = np.int(np.where(colRs=='y_f606w')[0])

    master_arr[:,0] = ref[:,xv]
    master_arr[:,1] = ref[:,yv]

    # Weights required for the 6D transformation function
    # We haven't assigned useful values to them.
    weights = np.ones(len(master_arr))

    # Creating an array that will hold the original and transformed values
    all_arr = np.zeros((len(all),2))

    # bt stands for before transform
    x_bt = np.int(np.where(colAs=='xcenter')[0])
    y_bt = np.int(np.where(colAs=='ycenter')[0])

    all_arr[:,0] = all[:,x_bt]
    all_arr[:,1] = all[:,y_bt]

    outName = dir + 'drcFiltTrans_' + targname + str(setnum)
    print('Transforming ',targname)
    # This takes the xy positions from the points you have matched
    # to the master, the xy positions of the corresponding points in
    # the master, weights, and the xy points in the match_arr frame
    # that will be transformed into the master_arr frame.
    new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1],
                                     master_arr[:,0],master_arr[:,1],weights,
                                     weights, all_arr[:,0],all_arr[:,1])

    outArr = np.hstack((all,new_all))  # attaching the transformed positions
    # to the original (F814W) catalog

    s0 = ' '
    header = s0.join(colAs)
    header += ' x_f606wTrans y_f606wTrans'

    np.savetxt(outName + '.dat',outArr,header=header)

    makePlot(targname,match_arr[:,0],match_arr[:,1],master_arr[:,0],
             master_arr[:,1],new_match[:,0],new_match[:,1],
             label_1='Original in F814W',label_2='Original in F606W',
             label_3='New in F814W 2 F606W',outname=outName+'_matchCheck')

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
