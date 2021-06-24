import numpy as np
import matplotlib.pyplot as plt
from hst_fun_2 import *
from scipy import stats

def magCorr(dataN):

    fileDir = '../Hannah_Data/'
    cor = hst_gc_z(fileDir+'wfc_F606W',dataN['xr1_f606w'],dataN['yr1_f606w'])

    cor2 = hst_gc_z(fileDir+'wfc_F814W',dataN['xr1_f814w'],dataN['yr1_f814w'])

    out = np.hstack((cor,cor2))

    return out

def outPutCorr(targname,dir='./'):

    fileN = np.genfromtxt(dir+targname+'_allMatchedZPTed_mDc.dat',names=True)

    file = np.genfromtxt(dir+targname+'_allMatchedZPTed_mDc.dat')

    colAs = np.array(fileN.dtype.names)

    s0=' '
    header = s0.join(colAs)

    corOut = np.zeros((len(fileN),2))
    cor = magCorr(fileN)

    # print(file.shape)
    # print(cor.shape)

    corOut[:,0] = cor[:,0]
    corOut[:,1] = cor[:,1]

    header += ' magCorr_f606w magCorr_f814w'

    out = np.hstack((file,corOut))

    np.savetxt(dir+'catWmagCorr.dat',out,header=header)

    return None


def plotMagCorrs(targname,filt,dir='./'):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    xstr = 'xt1'+fils
    ystr = 'yt1'+fils

    magStr = 'magCorr'+fils

    fileN = np.genfromtxt(dir+'catWmagCorr.dat',names=True)

    mean, x_edges, y_edges, _ = stats.binned_statistic_2d(fileN[xstr],fileN[ystr],fileN[magStr],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    cm = plt.cm.get_cmap('viridis')
    #
    bin_cent_y = (y_edges[1:] + y_edges[:-1])*0.5
    bin_cent_x = (x_edges[1:] + x_edges[:-1])*0.5
    #

    extent=(np.min(bin_cent_x),np.max(bin_cent_x),np.min(bin_cent_y),np.max(bin_cent_y))


    fig, ax = plt.subplots(figsize=(6,6))
    # #

    cax = ax.imshow(mean.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm)
    cbar = fig.colorbar(cax)

    title_str = targname + '_' + filt
    ax.set_title(title_str)

    plt.savefig(dir+targname+'_'+filt+'_binnedCorr.png',dpi=600,bbox_inches='tight')
    # plt.show()
    plt.close()

    return None
