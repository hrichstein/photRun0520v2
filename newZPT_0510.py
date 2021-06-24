import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def getZPT(targname,filt,dir='./',sigTol=2.5,stdTol=0.1):

    dir = './catRawMags20Aug/catDir_'+targname+'/'
    # drcDir = 'photUtils20Aug/catDir_'+targname+'/'
    file = np.genfromtxt(dir+targname+'_fd_'+filt+'_cut.dat',names=True)
    full = np.genfromtxt(dir+targname+'_flcDRCmatch_'+filt+'.dat',names=True)

    if filt=='F606W':
        fils='_f606w'
    elif filt=='F814W':
        fils='_f814w'

    magStr = 'magrD'+fils

    file_idx = file['stdF']<=stdTol
    file_g = file[file_idx]

    full_idx = full['stdF']<=stdTol
    full_g = full[full_idx]

    kG = True

    std_arr = file_g['stdF']
    std_full = full_g['stdF']

    drc_arr = file_g[magStr]
    drc_full = full_g[magStr]

    file_g = file_g[magStr]-file_g['magrF']
    full_g = full_g[magStr]-full_g['magrF']

    iter = 0
    while kG:
        print('Iteration',iter)
        # print('Length file_g',len(file_g))
        # print('Length std_arr',len(std_arr))
        #
        # print('Length full_g',len(full_g))
        # print('Length std_full',len(std_full))

        # w_avg_flc = sum(file_g * std_arr) / sum(std_arr)
        # w_avg_full = sum(full_g * std_full) / sum(std_full)

        w_avg_flc = np.nanmedian(file_g)
        w_avg_full = np.nanmedian(full_g)

        print('Weighted Average',w_avg_flc)

        temp_flc_diff = []
        temp_full_diff = []

        temp_drc_diff = []
        temp_drc_full = []

        temp_flc_err = []
        temp_full_err = []

        for cc,vv in enumerate(file_g):
            # if abs(vv/std_arr[cc]) <= abs(2.5 * w_avg_flc):
            if (w_avg_flc - (2.5*np.nanstd(file_g)) <= vv) and  (vv<= w_avg_flc + (2.5*np.nanstd(file_g))):
            # if (vv <= 2.5 * w_avg_flc) and (vv >= -2.5 * w_avg_flc):
                temp_flc_diff.append(vv)
                temp_flc_err.append(std_arr[cc])
                temp_drc_diff.append(drc_arr[cc])
        if len(temp_flc_diff)==0:
            temp_flc_diff = file_g
            temp_flc_err = std_arr
            temp_drc_diff = drc_arr

        for cc,vv in enumerate(full_g):
            # if abs(vv/std_full[cc]) <= abs(2.5 * w_avg_full):
            # if (vv <= 2.5 * w_avg_full) and (vv >= -2.5 * w_avg_full):
            if (w_avg_full - (2.5*np.nanstd(full_g)) <= vv) and (vv<= w_avg_full + (2.5*np.nanstd(full_g))):
                temp_full_diff.append(vv)
                temp_full_err.append(std_full[cc])
                temp_drc_full.append(drc_full[cc])
        if len(temp_full_diff)==0:
            temp_full_diff = full_g
            temp_full_err = std_full
            temp_drc_full = drc_full

        if (len(temp_flc_diff) == len(file_g)) and (len(temp_full_diff) == len(full_g)):
            kG = False

            full_diff = np.array(temp_full_diff)
            flc_diff = np.array(temp_flc_diff)

            drc_diff = np.array(temp_drc_diff)
            drc_full = np.array(temp_drc_full)

            std_arr = np.array(temp_flc_err)
            std_full = np.array(temp_full_err)

        else:
            # print('Before reassignment')
            #
            # print('Length temp_full_diff',len(temp_full_diff))
            # print('Length temp_full_err',len(temp_full_err))
            #
            # print('Length temp_flc_diff',len(temp_flc_diff))
            # print('Length temp_flc_err',len(temp_flc_err))
            kG = True
            file_g = np.array(temp_flc_diff)
            full_g = np.array(temp_full_diff)

            drc_arr = np.array(temp_drc_diff)
            drc_full = np.array(temp_drc_full)

            std_arr = np.array(temp_flc_err)
            std_full = np.array(temp_full_err)

            # print('After reassignment')
            # print('Length file_g',len(file_g))
            # print('Length std_arr',len(std_arr))
            #
            # print('Length full_g',len(full_g))
            # print('Length std_full',len(std_full))

            iter += 1

            if iter > 10:
                kG = False
                full_diff = np.array(temp_full_diff)
                flc_diff = np.array(temp_flc_diff)

                drc_diff = np.array(temp_drc_diff)
                drc_full = np.array(temp_drc_full)

                std_arr = np.array(temp_flc_err)
                std_full = np.array(temp_full_err)

    # print('Iteration',iter)
    # full_diff = np.array(full_diff)
    # flc_diff = np.array(flc_diff)
    #
    # std_arr = np.array(std_arr)
    # std_full = np.array(std_full)

    print('Length full',len(full_diff))
    print('Length cut',len(flc_diff))

    mag_corr = np.average(flc_diff,weights=1/(std_arr**2))
    med_corr = np.nanmedian(flc_diff)
    err_add = np.nanstd(flc_diff/ np.sqrt(len(flc_diff)))

    mag_full_corr = np.average(full_diff,weights=1/(std_full**2))
    med_full_corr = np.nanmedian(full_diff)
    err_full_add = np.nanstd(full_diff/ np.sqrt(len(full_diff)))

    # mag_corr = np.nanmedian(flc_diff[0])
    # err_add = np.nanstd(flc_diff[0] / np.sqrt(len(flc_diff[0])))

    # plotting part
    fig, ax = plt.subplots(figsize=(6,4))

    ax.scatter(drc_full,full_diff,s=20,color='silver',label='Full')
    ax.scatter(drc_diff,flc_diff,s=10,label='Cut')
    ax.hlines(mag_corr,18,26.5,label='Weight Avg Cut')
    ax.hlines(med_corr,18,26.5,label='Med Cut',color='cyan')
    ax.hlines(mag_full_corr,18,26.5,label='Weight Avg Full',color='hotpink')
    ax.hlines(med_full_corr,18,26.5,label='Med Full',color='lime')

    ax.set_xlim(18,26.5)
    title_str = targname + '_' + filt + ' {0}'.format(len(flc_diff))
    ax.set_ylabel('DRC - FLC (mag)')
    ax.set_xlabel('DRC (mag)')
    ax.set_title(title_str)

    ax.legend()

    # saveDir = './photUtils28Sep/catDir_'+targname+'/'
    saveDir = './photUtils21Oct/catDir_'+targname+'/'

    plt.savefig(saveDir+targname+'_'+filt+'ZPTline0510.png',dpi=600,bbox_inches='tight')

    plt.close()

    return mag_corr, med_corr,err_add, mag_full_corr, med_full_corr,err_full_add


def applyZPT(mag_corr, med_corr,err_add, mag_full_corr, med_full_corr,err_full_add, targname,filt,dir='./'):

    all = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    cat = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')

    mean = all['mean']
    std = all['stdev']

    new_mag = mean + mag_corr
    new_med = mean + med_corr
    new_err = np.sqrt ( std**2 + err_add**2)

    new_mag2 = mean + mag_full_corr
    new_med2 = mean + med_full_corr
    new_err2 = np.sqrt ( std**2 + err_full_add**2)

    newCol = np.zeros((len(all),6))
    newCol[:,0] = new_mag
    newCol[:,1] = new_med
    newCol[:,2] = new_err
    newCol[:,3] = new_mag2
    newCol[:,4] = new_med2
    newCol[:,5] = new_err2

    outArr = np.hstack((cat,newCol))

    colAs = np.array(all.dtype.names)
    s0=' '
    header = s0.join(colAs)
    header += ' magZPT_wa magZPT_med magZPTerr magZPTfull_wa magZPTfull_med magZPTerrfull'

    # saveDir = './photUtils28Sep/catDir_'+targname+'/'
    saveDir = './catRawMags20Aug/catDir_'+targname+'/'

    np.savetxt(saveDir+'magZPTedAll_'+filt+'.dat',outArr,header=header)

    return None


def doZPT(targname,filt,dir='./',sigTol=2.5,stdTol=0.1):

    mag_corr, med_corr,err_add, mag_full_corr, med_full_corr,err_full_add = getZPT(targname,filt,dir=dir,sigTol=sigTol,stdTol=stdTol)

    # applyZPT(mag_corr, med_corr,err_add, mag_full_corr, med_full_corr,err_full_add,targname,filt,dir=dir)

    applyZPT(mag_corr, med_corr,err_add, mag_full_corr, med_full_corr,err_full_add,targname,filt,dir='./catRawMags20Aug/catDir_'+targname+'/')
    return None
#
