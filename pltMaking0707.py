import numpy as np
import matplotlib.pyplot as plt
from getJdan import getJdan

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

def openFiles(targname,filt,dir='./',iter=1):

    jdanUse = getJdan(targname,filt)

    for dd in range(len(jdanUse)-1):
        file = np.genfromtxt(dir+jdanUse[dd+1]+"_"+targname+"_"+filt+"_t{0:d}.dat".format(iter),names=True)

        idx = np.argsort(file['magr'])[:50]

        file_cut = file[idx]

        match_arr = np.zeros((50,2))
        new_arr = np.zeros((50,2))

        match_arr[:,0] = file_cut['xo']
        match_arr[:,1] = file_cut['yo']

        new_arr[:,0] = file_cut['xt']
        new_arr[:,1] = file_cut['yt']

        makePlot(match_arr,new_arr,outname=dir+jdanUse[dd+1]+'_matchCheck.png')

    return None
