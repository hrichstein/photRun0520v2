import numpy as np
import matplotlib.pyplot as plt

def makeCMD(targname,dir='./'):

    drcDir = 'photUtils20Aug/catDir_'+targname+'/'
    drc = np.genfromtxt(drcDir+targname+'_filtMatchDRC_pU.dat',names=True)


    flc = np.genfromtxt(dir+targname+'_allMatchedZPTed_pu.dat',names=True)

    # d_idx = np.logical_and(drc['c_star_f606w']>=0.5, drc['c_star_f814w']>=0.5)
    # f_idx = np.logical_and(flc['c_star_f606w']>=0.5, flc['c_star_f814w']>=0.5)
    #
    # drc_g = drc[d_idx]
    # flc_g = flc[f_idx]

    drc_g = drc.copy()
    flc_g = flc.copy()

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)

    ax1.scatter(drc_g['magr_f606w']-drc_g['magr_f814w'],drc_g['magr_f606w'],s=5,label='DRC')
    ax1.scatter(flc_g['magZPT_f606w']-flc_g['magZPT_f814w'],flc_g['magZPT_f606w'],s=5,label='FLC')

    ax2.scatter(drc_g['magr_f606w']-drc_g['magr_f814w'],drc_g['magr_f814w'],s=5,label='DRC')
    ax2.scatter(flc_g['magZPT_f606w']-flc_g['magZPT_f814w'],flc_g['magZPT_f814w'],s=5,label='FLC')

    ax1.set_ylim(28,17.5)
    ax1.set_xlim(-1.5,1.5)
    ax1.set_xlabel('F606W-F814W')
    ax1.set_ylabel('F606W')

    ax2.set_xlabel('F606W-F814W')
    ax2.set_ylabel('F814W')

    ax1.set_title(targname)
    ax2.set_title(targname)

    ax1.legend()
    ax2.legend()

    plt.savefig(dir+targname+'_pu.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None
