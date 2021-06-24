""" Making CMDs to compare FLC/DRC CMDs and look at flags """

import numpy as np
import matplotlib.pyplot as plt


def make2CMD(targname,dir='./',saveDir='./'):

    file = np.genfromtxt(dir + targname + '_wErr.dat',names=True)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,7),
                                  sharey=True, sharex=True)

    # Full DRC CMD
    magr606 = -2.5 * np.log10(file['final_phot_f606w'])
    magr814 = -2.5 * np.log10(file['final_phot_f814w'])

    ax1.scatter(magr606-magr814,
                magr606,s=5,color='cornflowerblue')
    ax1.set_title('F606W')

    # Full DRC & FLC CMD

    ax2.scatter(magr606-magr814,
                magr814, s=5,
                color='rebeccapurple')

    # ax2.scatter(file['meanFLC_f606w']-file['meanFLC_f814w'],
    #             file['meanFLC_f606w'],s=5,label='FLC',color='darksalmon')
    ax2.set_title('F814W')

    # Flagged DRC & FLC CMD

    # mask = (file['six_4_flag_f606w']==1) \
    #     & (file['six_4_flag_f814w']==1)
    # cut = file[mask]

    # ax3.scatter(cut['final_phot_f606w']-cut['final_phot_f814w'],
    #             cut['final_phot_f606w'], s=5,label='DRC',color='cornflowerblue')
    #
    # ax3.scatter(cut['meanFLC_f606w']-cut['meanFLC_f814w'],
    #             cut['meanFLC_f606w'], s=5,label='FLC',color='darksalmon')
    # ax3.set_title('Cut on S/G DRC & FLC')

    ax1.set_ylim(2.5,-18)
    ax1.set_xlim(-1.55,1.55)
    ax1.set_ylabel('F606W')
    ax2.set_xlabel('F606W - F814W')

    ax2.set_ylabel('F814W')

    # ax1.legend()
    # ax2.legend()
    # ax3.legend()

    plt.subplots_adjust(hspace=0, wspace=0)

    plt.savefig(saveDir+targname+'_CMDsInstr.png',dpi=600,bbox_inches='tight')
    plt.close()

    return None
