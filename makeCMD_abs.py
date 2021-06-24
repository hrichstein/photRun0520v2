import numpy as np
import matplotlib.pyplot as plt

def makeCMDabs(targname,dir='./'):

    # drcDir = 'photUtils20Aug/catDir_'+targname+'/'
    # drc = np.genfromtxt(drcDir+targname+'_filtMatchDRC_pU.dat',names=True)
    flc = np.genfromtxt(dir+targname+'_absMag.dat',names=True)

    # drc_g = drc.copy()
    flc_g = flc.copy()

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)

    # ax1.scatter(drc_g['magr_f606w']-drc_g['magr_f814w'],drc_g['magr_f606w'],s=5,label='DRC')
    ax1.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f606w'],s=5,label='FLC')

    # ax2.scatter(drc_g['magr_f606w']-drc_g['magr_f814w'],drc_g['magr_f814w'],s=5,label='DRC')
    ax2.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f814w'],s=5,label='FLC')

    ax1.set_ylim(max(flc_g['magAbs_f814w']+0.5),min(flc_g['magAbs_f606w'])-0.5)
    ax1.set_xlim(-1.5,1.5)
    ax1.set_xlabel('F606W-F814W')
    ax1.set_ylabel('F606W')

    ax2.set_xlabel('F606W-F814W')
    ax2.set_ylabel('F814W')

    ax1.set_title(targname)
    ax2.set_title(targname)

    ax1.legend()
    ax2.legend()

    plt.savefig(dir+targname+'_cmdAbs.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None
