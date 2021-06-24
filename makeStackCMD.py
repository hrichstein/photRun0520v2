import numpy as np
import matplotlib.pyplot as plt

def f2mag_dirs(targname,date='20Aug',workDir='./'):

    return workDir+'catRawMags'+date+'/catDir_'+targname+'/'

def makeCMDabs(namefile):

    targname_arr = np.genfromtxt(namefile,dtype='str')

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)
    for tt,tname in enumerate(targname_arr):

        dir = f2mag_dirs(tname,date='20Aug',workDir='./')
        flc = np.genfromtxt(dir+tname+'_absMag.dat',names=True)

        flc_idx = np.logical_and(flc['magZPTerr_f606w']<0.025,flc['magZPTerr_f814w']<0.025)
        flc_g = flc[flc_idx]

        ax1.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f606w'],s=5)
        ax2.scatter(flc_g['magAbs_f606w']-flc_g['magAbs_f814w'],flc_g['magAbs_f814w'],s=5)

    ax1.set_ylim(12,-5)
    ax1.set_xlim(-1.5,1.5)
    ax1.set_xlabel('F606W-F814W')
    ax1.set_ylabel('F606W')

    ax2.set_xlabel('F606W-F814W')
    ax2.set_ylabel('F814W')

    plt.savefig('stacked_cmdAbs3.png',dpi=600,bbox_inches='tight')

    plt.close()


    return None


def main():

    makeCMDabs('targnamesDirections2.txt')


if __name__ == '__main__':
    main()

# targname_arr = np.genfromtxt('targnamesDirections2.txt',dtype='str')
