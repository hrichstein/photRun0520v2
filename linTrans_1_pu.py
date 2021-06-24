import numpy as np
from linear6d import *
from getJdan import getJdan
import matplotlib.pyplot as plt

def outDiths(targname,filt,jdanUse,dir='./',suffix='_ref.dat',iter=1):

    # jdanUse = getJdan(targname,filt)

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
        x_bt = np.int(np.where(colAs=='x_oc')[0])
        y_bt = np.int(np.where(colAs=='y_oc')[0])

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

            outArr = np.hstack((allCat,new_all))

            s0 = ' '
            header = s0.join(colAs)
            header += ' xt yt'

            np.savetxt(dir+outName,outArr,header=header)

        except RuntimeWarning:
            print('Not good enough.',targname,filt,dd)


    return None

def makePlot(match_arr,new_match,outname):

    match = match_arr
    new = new_match

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(match[:,0],match[:,1],s=30,label='Original')
    ax.scatter(new[:,0],new[:,1],s=15,label='New')

    ax.legend()
    ax.set_title(outname)

    # print(outname)
    # plt.show()
    plt.savefig(outname,dpi=600,bbox_inches='tight')
    #
    plt.close()

    return None

def openFiles(targname,filt,jdanUse,dir='./',iter=1):

    # jdanUse = getJdan(targname,filt)

    for dd in range(len(jdanUse)-1):
        file = np.genfromtxt(dir+jdanUse[dd+1]+"_"+targname+"_"+filt+"_t{0:d}.dat".format(iter),names=True)

        idx = np.argsort(file['magr'])[:50]

        file_cut = file[idx]

        match_arr = np.zeros((50,2))
        new_arr = np.zeros((50,2))

        match_arr[:,0] = file_cut['x_oc']
        match_arr[:,1] = file_cut['y_oc']

        new_arr[:,0] = file_cut['xt']
        new_arr[:,1] = file_cut['yt']

        makePlot(match_arr,new_arr,outname=dir+jdanUse[dd+1]+'_matchCheck.png')

    return None
#
