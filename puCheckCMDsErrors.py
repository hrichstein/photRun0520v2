import numpy as np
import matplotlib.pyplot as plt

def flcDRCpuCMD(targname):

    drc = np.genfromtxt('photUtils20Aug/catDir_'+targname+'/'+targname+'_drcCut.dat',names=True)

    flc = np.genfromtxt('catRawMags20Aug/catDir_'+targname+'/'+targname+'_allMatchedZPTed_pu.dat',names=True)


    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,7),sharey=True)

    # Plotting DRCs
    ax1.scatter(drc['magr_f606w']-drc['magr_f814w'],drc['magr_f606w'],label='DRC',s=5,color='black')
    ax2.scatter(drc['magr_f606w']-drc['magr_f814w'],drc['magr_f814w'],label='DRC',s=5,color='black')

    # Plotting FLCs
    ax1.scatter(flc['mean_f606w']-flc['mean_f814w'],flc['mean_f606w'],label='FLC',s=5,alpha=0.5,color='magenta')
    ax2.scatter(flc['mean_f606w']-flc['mean_f814w'],flc['mean_f814w'],label='FLC',s=5,alpha=0.5,color='magenta')

    ax1.set_ylabel('F606W')
    ax2.set_ylabel('F814W')

    ax1.set_xlabel('F606W-F814W')
    ax2.set_xlabel('F606W-F814W')

    ax1.set_ylim(28,17.5)

    ax1.set_xlim(-2,2)
    ax2.set_xlim(-2,2)

    ax1.legend()
    ax2.legend()

    plt.savefig('plots1210/'+targname+'_CMDcompDRCflc.png',dpi=600,bbox_inches='tight')

    return None

def flcErrPlots(targname):

    flc = np.genfromtxt('catRawMags20Aug/catDir_'+targname+'/'+targname+'_allMatchedZPTed_pu.dat',names=True)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,4))

    ax1.scatter(flc['mean_f606w'],flc['stdev_f606w'],s=5)
    ax2.scatter(flc['mean_f814w'],flc['stdev_f814w'],s=5)

    ax1.set_xlabel('F606W')
    ax2.set_xlabel('F814W')

    ax1.set_ylabel('STD F606W')
    ax2.set_ylabel('STD F814W')

    ax1.set_ylim(-0.01,0.25)
    ax2.set_ylim(-0.01,0.25)

    ax1.set_xlim(17.5,28)
    ax2.set_xlim(17.5,28)

    plt.savefig('plots1210/'+targname+'_errors.png',dpi=600,bbox_inches='tight')

    return None

def puCMDmultPoint(targname_arr,outname='TUCANA-II'):

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,7),sharey=True)

    for cc, targname in enumerate(targname_arr):
        drc = np.genfromtxt('photUtils20Aug/catDir_'+targname+'/'+targname+'_drcCut.dat',names=True)

        flc = np.genfromtxt('catRawMags20Aug/catDir_'+targname+'/'+targname+'_allMatchedZPTed_pu.dat',names=True)

        # Plotting DRCs
        ax1.scatter(drc['magr_f606w']-drc['magr_f814w'],drc['magr_f606w'],label='DRC',s=5,color='black')
        ax2.scatter(drc['magr_f606w']-drc['magr_f814w'],drc['magr_f814w'],label='DRC',s=5,color='black')

        # Plotting FLCs
        ax1.scatter(flc['mean_f606w']-flc['mean_f814w'],flc['mean_f606w'],label='FLC',s=5,alpha=0.5,color='magenta')
        ax2.scatter(flc['mean_f606w']-flc['mean_f814w'],flc['mean_f814w'],label='FLC',s=5,alpha=0.5,color='magenta')

    ax1.set_ylabel('F606W')
    ax2.set_ylabel('F814W')

    ax1.set_xlabel('F606W-F814W')
    ax2.set_xlabel('F606W-F814W')

    ax1.set_ylim(28,17.5)

    ax1.set_xlim(-2,2)
    ax2.set_xlim(-2,2)

    ax1.legend()
    ax2.legend()

    plt.savefig('plots1210/'+outname+'_fullCMDcompDRCflc.png',dpi=600,bbox_inches='tight')

    return None


# targname_arr = np.genfromtxt('targnamesPost.txt',dtype='str')
targname_arr = ['TUCANA-II-SE','TUCANA-II-NE']
# targname_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']

def main():
#
    for c1,targname in enumerate(targname_arr):
        # flcDRCpuCMD(targname)
        # flcErrPlots(targname)

        puCMDmultPoint(targname_arr,outname='TUCANA-II')

if __name__ == '__main__':
    main()
