import numpy as np
import matplotlib.pyplot as plt


def plotErr(targname,dir='./'):

    filename = targname + '_wErr.dat'
    fileN = np.genfromtxt(dir+filename,names=True)

    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,5),sharex=True,sharey=True)

    ax1.scatter(fileN['magr_f606w'],fileN['err_f606w'],s=5,
                color='mediumseagreen')
    ax2.scatter(fileN['magr_f814w'],fileN['err_f814w'],s=5,color='indianred')

    ax1.set_xlim(18,27)
    ax1.set_ylim(0,0.25)

    ax1.set_xlabel('F606W STMAG')
    ax2.set_xlabel('F814W STMAG')

    ax1.set_ylabel('F606W Err')
    ax2.set_ylabel('F814W Err')

    plt.savefig(dir+'magVerr_' + targname + '.png',dpi=600,
                bbox_inches='tight')

    plt.close()

    return None
