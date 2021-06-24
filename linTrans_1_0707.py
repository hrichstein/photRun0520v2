import numpy as np
from linear6d import *
from getJdan import getJdan
import matplotlib.pyplot as plt

def outDiths(targname,filt,dir='./',suffix='_ref.dat',iter=1):

    jdanUse = getJdan(targname,filt)

    file = np.genfromtxt(dir+'matched_w_MagsPos'+suffix,names=True)
    fileCat = np.genfromtxt(dir+'matched_w_MagsPos'+suffix)
    colNs = np.array(file.dtype.names)

    xt1 = np.int(np.where(colNs=='xt1')[0])
    yt1 = np.int(np.where(colNs=='yt1')[0])

    master_arr = np.zeros((len(file),2))

    master_arr[:,0] = fileCat[:,xt1]
    master_arr[:,1] = fileCat[:,yt1]

    match_arr = np.zeros((len(file),2))

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    for dd in range(len(jdanUse)-1):

        all = np.genfromtxt(dir+jdanUse[dd+1]+"_"+targname+"_"+filt+"_oc.dat",names=True)
        allCat = np.genfromtxt(dir+jdanUse[dd+1]+"_"+targname+"_"+filt+"_oc.dat")
        colAs = np.array(all.dtype.names)
        x_bt = np.int(np.where(colAs=='xo')[0])
        y_bt = np.int(np.where(colAs=='yo')[0])

        all_arr = np.zeros((len(all),2))

        all_arr[:,0] = allCat[:,x_bt]
        all_arr[:,1] = allCat[:,y_bt]

        outName = jdanUse[dd+1]+"_"+targname+"_"+filt+"_t{0:d}.dat".format(iter)

        str1 = 'xt'+str(dd+1)
        str2 = 'yt'+str(dd+1)

        xt = np.int(np.where(colNs==str1)[0])
        yt = np.int(np.where(colNs==str2)[0])

        match_arr[:,0] = fileCat[:,xt]
        match_arr[:,1] = fileCat[:,yt]

        try:
            new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

            # This would be a good place to have it make a plot
            # makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
            # new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Past FLC Trans Pos.',label_2='New in FLC2DRC',label_3='Original in DRC',outname=plotName)

            outArr = np.hstack((allCat,new_all))
            header = 'flags RA DEC xr yr flux c_star magr id xc yc xo yo xt yt'

            np.savetxt(dir+outName,outArr,header=header)

        except RuntimeWarning:
            print('Not good enough.',targname,filt,dd)


    return None

#
