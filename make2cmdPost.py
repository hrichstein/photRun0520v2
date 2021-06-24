import numpy as np
import matplotlib.pyplot as plt

def f2mag_dirs(targname,date='20Aug',workDir='./'):

    return workDir+'catRawMags'+date+'/catDir_'+targname+'/'

def makeNewCMD(sg_cut=0.5):

    targname_arr1 = ['HOROLOGIUM-I','PHOENIX-II','RETICULUM-II']
    targname_arr2 = ['HYDRA-II','PEGASUS-III','SAGITTARIUS-II',\
    'TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE',\
    'TUCANA-II-SE']

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)
    for tt,tname in enumerate(targname_arr1):

        dir = f2mag_dirs(tname,date='20Aug',workDir='./')
        flc = np.genfromtxt(dir+'classStarCat_F606Wmatch_'+tname+'.dat',names=True)

        flc_idx = np.logical_and(np.logical_and(flc['magZPT_f606w']>=19.5,flc['c_star_f606w']>=sg_cut),flc['c_star_f814w']>=sg_cut)
        flc_g = flc[flc_idx]

        ax1.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f606w'],s=5)

    for tt,tname in enumerate(targname_arr2):

        dir = f2mag_dirs(tname,date='20Aug',workDir='./')
        flc = np.genfromtxt(dir+'classStarCat_F606Wmatch_'+tname+'.dat',names=True)

        flc_idx = np.logical_and(np.logical_and(flc['magZPT_f606w']>=19.5,flc['c_star_f606w']>=sg_cut),flc['c_star_f814w']>=sg_cut)
        flc_g = flc[flc_idx]

        ax2.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f606w'],s=5)
        # ax2.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f814w'],s=5)

    ax1.set_ylim(12,-5)
    ax1.set_xlim(-1.5,1.5)
    ax1.set_xlabel('F606W-F814W')
    ax1.set_ylabel('F606W')
    ax1.set_title('LMC Associated')

    ax2.set_xlabel('F606W-F814W')
    ax2.set_ylabel('F606W')
    ax2.set_title('MW Long-Term')

    # plt.savefig('catRawMags20Aug/stacked_cmdAbs_{0}.png'.format(str(sg_cut)),dpi=600,bbox_inches='tight')
    plt.savefig('catRawMags20Aug/poster_cmdSplit_{0}.png'.format(str(sg_cut)),dpi=600,bbox_inches='tight')

    plt.close()

    return None

def main():

    # makeNewCMD('targnamesDirections2.txt',sg_cut=0.95)
    makeNewCMD(sg_cut=0.95)


if __name__ == '__main__':
    main()
