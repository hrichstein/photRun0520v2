import numpy as np
import matplotlib.pyplot as plt
from hst_fun_2 import *

from getZPT0707 import *
from mpl_toolkits.axes_grid1 import make_axes_locatable


def magCorr(dataN,filt):

    fileDir = '../Hannah_Data/'
    cor = hst_gc_z(fileDir+'wfc_'+filt,dataN['xR'],dataN['yR'])


    return cor

def magCorrTry(targname,filt,dir='./'):

    mag_corr, err_add = getZPT(targname,filt,dir=dir,sigTol=2.5,stdTol=0.05)

    fileN = np.genfromtxt(dir+'flcDRCmatch_'+filt+'_mDc.dat',names=True)

    file = np.genfromtxt(dir+'flcDRCmatch_'+filt+'_mDc.dat')

    colNs = np.array(fileN.dtype.names)

    mean_idx = np.int(np.where(colNs=='magrF')[0])
    mean = file[:,mean_idx]

    zpt_mag = np.zeros((len(mean),1))
    zpt_mag[:,0] = mean + mag_corr

    cor = magCorr(fileN,filt)
    cor_mag = zpt_mag + cor

    cor2 = zpt_mag + cor - mag_corr

    fileOut = np.hstack((file,zpt_mag,cor,cor_mag,cor2))

    s0 = ' '
    header = s0.join(colNs)
    header += ' magZPT dCorr magDcorr magDcorrNZPT'
    #
    # form = '%1.5f %1.5f %1.4f %1.5f %d %1.5f %1.5f %d %1.5f %1.5f %1.4f %1.4f %1.4f %1.4f %1.4f '

    np.savetxt(dir+'magCorrFile_'+filt+'.dat',fileOut,header=header)

    return None

def plotFunc(targname,filt,dir='./',sigTol=2.5,stdTol=0.05):

    file = np.genfromtxt(dir+'magCorrFile_'+filt+'.dat',names=True)

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

    # magDiff = file_g['magrD']-file_g['magZPT']
    magDiff = file_g['magrD']-file_g['magrF']
    use = np.logical_and(magDiff>=min_diff,magDiff<=max_diff)

    meanOrig, x_edgesOrig, y_edgesOrig, _ = stats.binned_statistic_2d(file_g['xD'][use],file_g['yD'][use],magDiff[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    meanCorrVal, x_edgesCorrVal, y_edgesCorrVal, _ = stats.binned_statistic_2d(file_g['xR'][use],file_g['yR'][use],file_g['dCorr'][use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    magDiff2 = file_g['magrD']-file_g['magDcorrNZPT']

    meanMagCorr, x_edgesMagCorr, y_edgesMagCorr, _ = stats.binned_statistic_2d(file_g['xD'][use],file_g['yD'][use],magDiff2[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    magDiff3 = file_g['magrD']-file_g['magDcorr']
    meanMagCorrZ, x_edgesMagCorrZ, y_edgesMagCorrZ, _ = stats.binned_statistic_2d(file_g['xD'][use],file_g['yD'][use],magDiff3[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    magDiff4 = file_g['magrD']-file_g['magZPT']
    meanZDiff, x_edgesMagCorrZ, y_edgesMagCorrZ, _ = stats.binned_statistic_2d(file_g['xD'][use],file_g['yD'][use],magDiff4[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])


    min_diff2 = np.min(meanZDiff)
    max_diff2 = np.max(meanZDiff)

    bin_cent_y = (y_edgesOrig[1:] + y_edgesOrig[:-1])*0.5
    bin_cent_x = (x_edgesOrig[1:] + x_edgesOrig[:-1])*0.5
    extent=(np.min(bin_cent_x),np.max(bin_cent_x),np.min(bin_cent_y),np.max(bin_cent_y))

    cm = plt.cm.get_cmap('viridis')
    #
    fig, ((ax1,ax2), (ax3,ax4) )= plt.subplots(2,2,figsize=(12,10))

    cax1 = ax1.imshow(meanOrig.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm)
    # cbar1 = fig.colorbar(cax1)
    # ax1.set_title(targname+'_'+filt+'_'+'DRC - ZPT_FLC')
    ax1.set_title(targname+'_'+filt+'_'+'DRC - FLC')
    ax1.set_xlabel('X Pixel (DRC)')
    ax1.set_ylabel('Y Pixel (DRC)')

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cax1, cax=cax, orientation='vertical')
    # cax2 = ax2.imshow(meanCorrVal.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm)
    # cbar2 = fig.colorbar(cax2)
    # ax2.set_title(targname+'_'+filt+'_'+'Correction')
    # ax2.set_xlabel('X Pixel (DRC)')
    # ax2.set_ylabel('Y Pixel (DRC)')

    cax2 = ax2.imshow(meanMagCorr.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm)
    # cbar2 = fig.colorbar(cax2)
    ax2.set_title('DRC - DC_FLC')
    ax2.set_xlabel('X Pixel (DRC)')
    ax2.set_ylabel('Y Pixel (DRC)')

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cax2, cax=cax, orientation='vertical')

    cax3 = ax3.imshow(meanZDiff.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm,vmin=min_diff2,vmax=max_diff2)
    # cbar3 = fig.colorbar(cax3)
    ax3.set_title('DRC - ZPT_FLC')
    ax3.set_xlabel('X Pixel (DRC)')
    ax3.set_ylabel('Y Pixel (DRC)')

    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cax3, cax=cax, orientation='vertical')

    cax4 = ax4.imshow(meanMagCorrZ.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm,vmin=min_diff2,vmax=max_diff2)
    # cbar4 = fig.colorbar(cax4)
    ax4.set_title('DRC - ZPT_DC_FLC')
    ax4.set_xlabel('X Pixel (DRC)')
    ax4.set_ylabel('Y Pixel (DRC)')

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cax4, cax=cax, orientation='vertical')


    plt.savefig(dir+targname+'_'+filt+'ZPTbinned_comp_Add.png',dpi=600,bbox_inches='tight')
    # plt.show()
    plt.close()

    # fig,ax = plt.subplots(figsize=(6,6))
    #
    # cax1 = ax.imshow(meanCorrVal.T,extent=extent,origin='lower',interpolation='nearest',cmap=cm)
    # # cbar1 = fig.colorbar(cax1)
    # # ax1.set_title(targname+'_'+filt+'_'+'DRC - ZPT_FLC')
    # ax.set_title(targname+'_'+filt+'_'+'Correction')
    # ax.set_xlabel('X Pixel (DRC)')
    # ax.set_ylabel('Y Pixel (DRC)')
    #
    # # divider = make_axes_locatable(ax1)
    # # cax = divider.append_axes('right', size='5%', pad=0.05)
    # # fig.colorbar(cax1, cax=cax, orientation='vertical')
    # cbar1 = fig.colorbar(cax1)
    #
    # plt.savefig(dir+targname+'_'+filt+'ZPTbinned_corrVal.png',dpi=600,bbox_inches='tight')
    # # plt.show()
    # plt.close()


    return None

def doAll(targname,filt,dir='./'):

    magCorrTry(targname,filt,dir=dir)
    plotFunc(targname,filt,dir=dir)

    return None


    #
