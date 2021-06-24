import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def plotZPTdiff(targname,filt,dir='./',sigTol=2.5,stdTol=0.05):


    # file = np.genfromtxt(dir+'flcDRCmatch_'+filt+'.dat',names=True)
    file = np.genfromtxt(dir+'flcDRCmatch_'+filt+'_mDc.dat',names=True)
    drc = np.genfromtxt(dir+'drc_useful_'+targname+'.dat')

    lenD = len(drc)

    kG = True # keep going

    # Loop to try to make sure I'm basing the magnitude correction off of a decent amount of stars.
    while kG:

        file_g = file[file['stdF']<=stdTol]
        # file_g = file[file_idx]

        flc_diff = stats.sigmaclip(file_g['magrD']-file_g['magrF'],sigTol,sigTol)

        if len(flc_diff[0]) >= (0.1*lenD):
            kG = False

        else:
            stdTol += 0.05
            sigTol += 0.25

        if (stdTol >= 0.5) or (sigTol >= 3.5):
            print('Not very good stats. Need more stars.')
            kG = False

    min_diff = np.min(flc_diff[0])
    max_diff = np.max(flc_diff[0])

    magDiff = file_g['magrD']-file_g['magrF']
    use = np.logical_and(magDiff>=min_diff,magDiff<=max_diff)

    mean, x_edges, y_edges, _ = stats.binned_statistic_2d(file_g['xD'][use],file_g['yD'][use],magDiff[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])


    # plotting part
    cm = plt.cm.get_cmap('viridis')
    #
    bin_cent_y = (y_edges[1:] + y_edges[:-1])*0.5
    bin_cent_x = (x_edges[1:] + x_edges[:-1])*0.5
    #

    extent=(np.min(bin_cent_x),np.max(bin_cent_x),np.min(bin_cent_y),np.max(bin_cent_y))
    # xv, yv = np.meshgrid(bin_cent_x,bin_cent_y)
    #

    fig, ax = plt.subplots(figsize=(6,6))
    # #

    cax = ax.imshow(mean.T,extent=extent,origin='lower',interpolation='nearest',vmin=min_diff,vmax=max_diff,cmap=cm)
    cbar = fig.colorbar(cax)
    # # sc = ax.scatter(xv.flatten(),yv.flatten(),c=mean.flatten(),vmin=min_diff,vmax=max_diff,s=20,cmap=cm)
    #
    # ax.pcolormesh(xv,yv,mean.T,snap=True)
    # # # sc = ax.scatter(file_g['xD'][use],file_g['yD'][use],c=magDiff[use],vmin=min_diff,vmax=max_diff,s=20,cmap=cm)
    # # #
    # plt.colorbar(sc)
    title_str = targname + '_' + filt + ' {0}'.format(len(flc_diff[0]))
    ax.set_title(title_str)
    # #
    # plt.savefig(dir+targname+'_'+filt+'ZPTbinned.png',dpi=600,bbox_inches='tight')
    plt.savefig(dir+targname+'_'+filt+'ZPTbinned_mDc.png',dpi=600,bbox_inches='tight')
    # plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(6,6))
    sc = ax.scatter(file_g['xD'][use],file_g['yD'][use],c=magDiff[use],vmin=min_diff,vmax=max_diff,s=15,cmap=cm)

    ax.set_xlim(np.min(bin_cent_x),np.max(bin_cent_x))
    ax.set_ylim(np.min(bin_cent_y),np.max(bin_cent_y))
    plt.colorbar(sc)
    title_str = targname + '_' + filt + ' {0}'.format(len(flc_diff[0]))
    ax.set_title(title_str)
    # #
    # plt.savefig(dir+targname+'_'+filt+'ZPTscatter.png',dpi=600,bbox_inches='tight')
    plt.savefig(dir+targname+'_'+filt+'ZPTscatter_mDc.png',dpi=600,bbox_inches='tight')
    # plt.show()
    plt.close()


    return None
