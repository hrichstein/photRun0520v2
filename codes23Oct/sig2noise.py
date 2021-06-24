import numpy as np
import matplotlib.pyplot as plt


def s2n(targname,dir='./'):

    filename = targname + '_wErr.dat'

    fileN = np.genfromtxt(dir+filename,names=True)
    file = np.genfromtxt(dir+filename)

    colFs = np.array(fileN.dtype.names)

    err_f606w = np.int(np.where(colFs=='err_f606w')[0])
    err_f814w = np.int(np.where(colFs=='err_f814w')[0])

    newCols = np.zeros((len(file),2))
    newCols[:,0] = 1.0875 * (1/file[:,err_f606w])
    newCols[:,1] = 1.0875 * (1/file[:,err_f814w])

    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,5),sharex=True,sharey=True)

    ax1.scatter(fileN['magr_f606w'],newCols[:,0],s=5,color='seagreen')
    ax2.scatter(fileN['magr_f814w'],newCols[:,1],s=5,color='indianred')

    ax1.set_title('{0} F606W'.format(targname))
    ax2.set_title('{0} F814W'.format(targname))

    ax1.set_ylabel('S/N')

    ax1.set_xlabel('F606W STMAG')
    ax2.set_xlabel('F814W STMAG')

    ax1.set_xlim(18,28)
    ax1.set_ylim(100,0)

    # plt.subplots_adjust(wspace=0,hspace=0)
    plt.savefig(dir+'sig2noise_' + targname + '.png',dpi=600,
                bbox_inches='tight')

    plt.close()

    return None

#
