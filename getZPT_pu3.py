import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def getZPT(targname,filt,dir='./',sigTol=2.5,stdTol=0.05):


    # file = np.genfromtxt(dir+'flcDRCmatch_'+filt+'.dat',names=True)
    # drc = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat')
    #
    # lenD = len(drc)


    file = np.genfromtxt(dir+targname+'_fd_'+filt+'_cut.dat',names=True)
    full = np.genfromtxt(dir+targname+'_flcDRCmatch_'+filt+'.dat',names=True)

    if filt=='F606W':
        fils='_f606w'
    elif filt=='F814W':
        fils='_f814w'

    magStr = 'magrD'+fils

    kG = True # keep going

    # Loop to try to make sure I'm basing the magnitude correction off of a decent amount of stars.
    while kG:

        file_g = file[file['stdF']<=stdTol]
        # file_g = file[file_idx]

        flc_diff = stats.sigmaclip(file_g[magStr]-file_g['magrF'],sigTol,sigTol)

        if len(flc_diff[0]) >= (0.1*len(full)):
            kG = False

        else:
            stdTol += 0.05
            sigTol += 0.25

        if (stdTol >= 0.5) or (sigTol >= 3.5):
            print('Not very good stats. Need more stars.')
            kG = False

    mag_corr = np.nanmedian(flc_diff[0])
    err_add = np.nanstd(flc_diff[0] / np.sqrt(len(flc_diff[0])))

    # plotting part

    fig, ax = plt.subplots(figsize=(6,4))

    ax.scatter(full[magStr],full[magStr]-full['magrF'],s=15,color='silver')
    ax.scatter(file_g[magStr],file_g[magStr]-file_g['magrF'],s=20)
    ax.hlines(mag_corr,17.5,28)

    ax.set_xlim(17.5,28)
    title_str = targname + '_' + filt + ' {0}'.format(len(flc_diff[0]))
    ax.set_title(title_str)

    saveDir = './photUtils28Sep/catDir_'+targname+'/'

    plt.savefig(saveDir+targname+'_'+filt+'ZPTline3.png',dpi=600,bbox_inches='tight')

    plt.close()

    return mag_corr, err_add


def applyZPT(mag_corr,err_add,targname,filt,dir='./'):

    all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    cat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')

    mean = all['mean']
    std = all['stdev']

    new_mag = mean + mag_corr
    new_err = np.sqrt ( std**2 + err_add**2)

    newCol = np.zeros((len(all),2))
    newCol[:,0] = new_mag
    newCol[:,1] = new_err

    outArr = np.hstack((cat,newCol))

    colAs = np.array(all.dtype.names)
    s0=' '
    header = s0.join(colAs)
    header += ' magZPT magZPTerr'

    saveDir = './photUtils28Sep/catDir_'+targname+'/'

    np.savetxt(saveDir+'magZPTedAll_'+filt+'.dat',outArr,header=header)

    return None


def doZPT(targname,filt,dir='./',sigTol=2.5,stdTol=0.05):

    mag_corr,err_add = getZPT(targname,filt,dir=dir,sigTol=2.5,stdTol=stdTol)

    applyZPT(mag_corr,err_add,targname,filt,dir=dir)


    return None
#
