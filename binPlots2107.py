import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

file = np.genfromtxt('./catRawMags1305/catDir_HOROLOGIUM-I_pastRuns/matchedFLCpsf2506_all.dat',names=True)

sigTol = 3.5

magDiff_aperDRC_606 = file['magDRC_f606w']-file['m606cAPER']
magDiff_aperDRC_814 = file['magDRC_f814w']-file['m814cAPER']

magDiff_psfDRC_606 = file['magDRC_f606w']-file['m606cPSF']
magDiff_psfDRC_814 = file['magDRC_f814w']-file['m814cPSF']

magDiff_flcDRC_606 = file['magDRC_f606w']-file['magZPT_f606w']
magDiff_flcDRC_814 = file['magDRC_f814w']-file['magZPT_f814w']

magDiff_flcAPER_606 = file['m606cAPER']-file['magZPT_f606w']
magDiff_flcAPER_814 = file['m814cAPER']-file['magZPT_f814w']

magDiff_flcPSF_606 = file['m606cPSF']-file['magZPT_f606w']
magDiff_flcPSF_814 = file['m814cPSF']-file['magZPT_f814w']

magDiff_aperPSF_606 = file['m606cAPER']-file['m606cPSF']
magDiff_aperPSF_814 = file['m814cAPER']-file['m814cPSF']

x_606 = file['xDRC_mat_f606w']
y_606 = file['yDRC_mat_f606w']

x_814 = file['xDRC_mat_f814w']
y_814 = file['yDRC_mat_f814w']

def makeBinPlot(x,y,diff_array,filt,outname,sigTol=3.5):

    cut_arr = stats.sigmaclip(diff_array,sigTol,sigTol)

    min_diff = np.min(cut_arr[0])
    max_diff = np.max(cut_arr[0])

    use = np.logical_and(diff_array>=min_diff,diff_array<=max_diff)

    mean, x_edges, y_edges, _ = stats.binned_statistic_2d(x[use],y[use],diff_array[use],statistic='mean',bins=20,range=[[0,4096],[0,4096]])

    cm = plt.cm.get_cmap('viridis')

    bin_cent_y = (y_edges[1:] + y_edges[:-1])*0.5
    bin_cent_x = (x_edges[1:] + x_edges[:-1])*0.5
    #

    extent=(np.min(bin_cent_x),np.max(bin_cent_x),np.min(bin_cent_y),np.max(bin_cent_y))

    fig, ax = plt.subplots(figsize=(6,6))

    cax = ax.imshow(mean.T,extent=extent,origin='lower',interpolation='nearest',vmin=min_diff,vmax=max_diff,cmap=cm)
    cbar = fig.colorbar(cax)

    title_str =  outname + filt
    ax.set_title(title_str)

    plt.savefig('./plots2107/' +outname+filt+'_bin.png',dpi=600,bbox_inches='tight')

    return None


# makeBinPlot(x_606,y_606,magDiff_aperDRC_606,'F606W','DRC-APER_',sigTol=3.5)
#
# makeBinPlot(x_606,y_606,magDiff_psfDRC_606,'F606W','DRC-PSF_',sigTol=3.5)
#
# makeBinPlot(x_606,y_606,magDiff_flcDRC_606,'F606W','DRC-FLC_',sigTol=3.5)
#
# makeBinPlot(x_606,y_606,magDiff_flcAPER_606,'F606W','APER-FLC_',sigTol=3.5)
#
# makeBinPlot(x_606,y_606,magDiff_flcPSF_606,'F606W','PSF-FLC_',sigTol=3.5)
#
# #F814W
# makeBinPlot(x_814,y_814,magDiff_aperDRC_814,'F814W','DRC-APER_',sigTol=3.5)
#
# makeBinPlot(x_814,y_814,magDiff_psfDRC_814,'F814W','DRC-PSF_',sigTol=3.5)
#
# makeBinPlot(x_814,y_814,magDiff_flcDRC_814,'F814W','DRC-FLC_',sigTol=3.5)
#
# makeBinPlot(x_814,y_814,magDiff_flcAPER_814,'F814W','APER-FLC_',sigTol=3.5)
#
# makeBinPlot(x_814,y_814,magDiff_flcPSF_814,'F814W','PSF-FLC_',sigTol=3.5)
makeBinPlot(x_814,y_814,magDiff_aperPSF_814,'F814W','APER-PSF_',sigTol=3.5)

makeBinPlot(x_606,y_606,magDiff_aperPSF_606,'F606W','APER-PSF_',sigTol=3.5)


#
