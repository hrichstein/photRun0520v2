""" Making CMDs to compare FLC/DRC CMDs and look at flags """

import numpy as np
import matplotlib.pyplot as plt


def make2CMD(targname,dir='./'):

    file = np.genfromtxt(dir + targname + '_wErr.dat',names=True)

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,7),
                                      sharey=True, sharex=True)

    # Full DRC CMD

    ax1.scatter(file['magr_f606w']-file['magr_f814w'],
                file['magr_f606w'],s=5,label='DRC',color='cornflowerblue')
    ax1.set_title('DRC Only')

    # Full DRC & FLC CMD

    ax2.scatter(file['magr_f606w']-file['magr_f814w'],
                file['magr_f606w'], s=5,label='DRC',color='cornflowerblue')

    ax2.scatter(file['meanFLC_f606w']-file['meanFLC_f814w'],
                file['meanFLC_f606w'],s=5,label='FLC',color='darksalmon')
    ax2.set_title('DRC & FLC')

    # Flagged DRC & FLC CMD

    mask = np.logical_or(file['six_4_flag_f606w']==1,
                         file['six_4_flag_f814w']==1)
    cut = file[mask]

    ax3.scatter(cut['magr_f606w']-cut['magr_f814w'],
                cut['magr_f606w'], s=5,label='DRC',color='cornflowerblue')

    ax3.scatter(cut['meanFLC_f606w']-cut['meanFLC_f814w'],
                cut['meanFLC_f606w'], s=5,label='FLC',color='darksalmon')
    ax3.set_title('Cut on S/G DRC & FLC')

    ax1.set_ylim(27.5,18.5)
    ax1.set_xlim(-1.55,1.55)
    ax1.set_ylabel('F606W')
    ax2.set_xlabel('F606W - F814W')

    ax1.legend()
    ax2.legend()
    ax3.legend()

    plt.subplots_adjust(hspace=0, wspace=0)

    plt.savefig(dir+targname+'_CMDs.png',dpi=600,bbox_inches='tight')
    plt.close()

    return None
