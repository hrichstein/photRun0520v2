""" Making CMDs to compare photometry apertures/techniques """

import numpy as np
import matplotlib.pyplot as plt


def main():
    targname_arr = np.genfromtxt('./targnamesPost.txt', dtype=str)
    # targname='HOROLOGIUM-I'
    # save_dir = './catMatchFLCdrc18Oct/catDir_' + targname + '/'
    # make3CMD(targname,saveDir=save_dir)
    for c1, targname in enumerate(targname_arr):
        save_dir = './catMatchFLCdrc18Oct/catDir_' + targname + '/'
        make3CMD(targname,saveDir=save_dir)


def make3CMD(targname,saveDir='./'):

    # SE and Pu3 file will be the same
    # Naming the files I need
    se_dir = './catMatchFLCdrc18Oct/catDir_' + targname + '/'
    se_file = se_dir + targname + '_fullSEpu.dat'
    # pu4_file = se_file
    # se_dir = './catMatchFLCdrc18Oct/seDRCs/'
    # se_file = se_dir + targname + '_sfErr.dat'
    # # magRaw_v, magRaw_i
    #
    pu3_dir = './photUtils21Oct/catDir_' + targname + '/'
    pu_3pix = pu3_dir + targname + '_filtMatchDRC_pU.dat'
    # # mean_f606w, mean_f814w, stdev_f606w, stdev_f814w

    pu4dir = './catRawMags20Aug/catDir_' + targname + '/'
    pu4pixE_f = pu4dir + targname + '_allMatchedZPTed_pu.dat'

    pu3dir = './photUtils21Oct/catDir_' + targname + '/'
    pu3pixE_f = pu3dir + targname + '_allMatchedZPTed_pu.dat'

    # mean_f606w, mean_f814w, stdev_f606w, stdev_f814w

    # Importing them
    se_3pix = np.genfromtxt(se_file,names=True)
    pu_3pix = np.genfromtxt(pu_3pix,names=True)
    pu_4pix = np.genfromtxt(se_file,names=True)
    pu4pix_E = np.genfromtxt(pu4pixE_f,names=True)
    pu3pix_E = np.genfromtxt(pu3pixE_f,names=True)

    # Plotting time
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,9))

    # 3 Pixel: SE vs PU
    ax1.scatter(se_3pix['magRaw_v']
                - se_3pix['magRaw_i'],se_3pix['magRaw_v'],s=5,
                label='SE 3 Pix',color='blue')
    ax1.scatter(pu_3pix['magr_f606w']
                -pu_3pix['magr_f814w'],pu_3pix['magr_f814w'],s=5,
                label='PU 3 Pix',color='magenta')
    ax1.set_title('3 Pixel: SE vs PU')
    ax1.set_ylim(26.5,18)
    ax1.set_xlim(-1.5,1.25)

    # PU: 3 Pixel vs 4 Pixel
    ax2.scatter(pu_4pix['magr_f606w']
                - pu_4pix['magr_f814w'],pu_4pix['magr_f606w'],s=5,
                label='PU 4 Pix',color='orange')
    ax2.scatter(pu_3pix['magr_f606w']
                - pu_3pix['magr_f814w'],pu_3pix['magr_f814w'],s=5,
                label='PU 3 Pix',color='magenta')
    ax2.set_title('PU: 3 Pixel vs 4 Pixel')
    ax2.set_ylim(26.5,18)
    ax2.set_xlim(-1.5,1.25)

    # SE 3 Pix vs PU 4 Pixel
    ax3.scatter(se_3pix['magRaw_v']
                - se_3pix['magRaw_i'],se_3pix['magRaw_v'],s=5,
                label='SE 3 Pix',color='blue')
    ax3.scatter(pu_4pix['magr_f606w']
                - pu_4pix['magr_f814w'],pu_4pix['magr_f606w'],s=5,
                label='PU 4 Pix',color='orange')
    ax3.set_title('SE 3 Pix vs PU 4 Pixel')
    ax3.set_ylim(26.5,18)
    ax3.set_xlim(-1.5,1.25)

    # PU: 3 Pix vs 4 Pix STD
    # ax4.scatter(pu4pix_E['mean_f606w'],pu4pix_E['stdev_f606w'],s=5,
    #             label='4 Pix PU',color='black')
    # ax4.scatter(pu3pix_E['mean_f606w'],pu3pix_E['stdev_f606w'],s=5,
    #             label='3 Pix PU',alpha=0.5,color='pink')
    ax4.set_xlim(18,26.5)
    ax4.set_ylim(-1e-4,0.3)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    # ax4.legend()

    plt.savefig(saveDir+'compCMDandSTD.png',dpi=600,bbox_inches='tight')
    plt.close()
    # print('Mean 3 Pix STD:',np.nanmean(pu3pix_E['stdev_f606w']))
    # print('Median 3 Pix STD:',np.nanmedian(pu3pix_E['stdev_f606w']))
    # print('Mean 4 Pix STD:',np.nanmean(pu4pix_E['stdev_f606w']))
    # print('Median 4 Pix STD:',np.nanmedian(pu4pix_E['stdev_f606w']))

    return None


if __name__ == '__main__':
    main()

#
