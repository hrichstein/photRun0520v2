import numpy as np
import matplotlib.pyplot as plt

# def f2mag_dirs(targname,date='20Aug',workDir='./'):
#
#     return workDir+'catRawMags'+date+'/catDir_'+targname+'/'
catDir = 'drcPhot10Nov/'

# plotDir = catDir + 'paulCMDs/'

tri_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']
tuc_arr = ['TUCANA-II-NE','TUCANA-II-SE']
urs_arr = ['URSA-MAJOR-II-WEST','URSA-MAJOR-II-EAST']
tu4_arr = ['TUCANA-IV-SOUTH','TUCANA-IV-NORTH']
tu3_arr = ['TUCANA-III-WEST','TUCANA-III-EAST']
seg_arr = ['SEGUE-1-WEST','SEGUE-1-EAST']
hor_arr = ['HOROLOGIUM-II-WEST','HOROLOGIUM-II-EAST']
boo_arr = ['BOOTES-II-SOUTH','BOOTES-II-NORTH']


def concat(name_arr,source,catDir='../drcPhot10Nov/'):

    for ff in range(len(name_arr)):

        if ff==0:
            f_out = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                  + name_arr[ff] + '_fullCat.dat')
            fhead = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                  + name_arr[ff] + '_fullCat.dat',
                                  names=True)
        else:
            f_temp = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                   + name_arr[ff] + '_fullCat.dat')
            f_out = np.vstack((f_out,f_temp))

    colNs = np.array(fhead.dtype.names)
    header = ' '.join(colNs)

    np.savetxt(catDir+ 'catDir_' + name_arr[ff] + '/' + source
               + "_fullCat_ALL.dat",f_out,header=header,fmt='%1.5f')

    makeCMD(catDir+ 'catDir_' + name_arr[ff] + '/',source)

    return None


def makeCMD(dir,source):

    file = np.genfromtxt(dir + source + "_fullCat_ALL.dat",names=True)

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,7),
                                      sharey=True, sharex=True)

    # Full DRC CMD

    ax1.scatter(file['magr_f606w']-file['magr_f814w'],
                file['magr_f606w'],s=5,label='DRC',color='cornflowerblue')
    ax1.set_title('DRC Only')

    # Full DRC & FLC CMD

    ax2.scatter(file['magr_f606w']-file['magr_f814w'],
                file['magr_f606w'], s=5,label='DRC',color='cornflowerblue')

    ax2.scatter(file['mean_f606w']-file['mean_f814w'],
                file['mean_f606w'],s=5,label='FLC',color='darksalmon')
    ax2.set_title('DRC & FLC')

    # Flagged DRC & FLC CMD

    mask = (file['six_4_flag_f606w']==1) \
        & (file['six_4_flag_f814w']==1)
    cut = file[mask]

    ax3.scatter(cut['magr_f606w']-cut['magr_f814w'],
                cut['magr_f606w'], s=5,label='DRC',color='cornflowerblue')

    ax3.scatter(cut['mean_f606w']-cut['mean_f814w'],
                cut['mean_f606w'], s=5,label='FLC',color='darksalmon')
    ax3.set_title('Cut on S/G DRC & FLC')

    ax1.set_ylim(27.5,18.5)
    ax1.set_xlim(-1.55,1.55)
    ax1.set_ylabel('F606W')
    ax2.set_xlabel('F606W - F814W')

    ax1.legend()
    ax2.legend()
    ax3.legend()

    plt.subplots_adjust(hspace=0, wspace=0)

    plt.savefig(dir + source + '_allCMDs.png',dpi=600,
                bbox_inches='tight')

    plt.close()

    return None


concat(tri_arr,'TRIANGULUM-II')
concat(tuc_arr,'TUCANA-II')
concat(urs_arr,'URSA-MAJOR-II')
concat(tu4_arr,'TUCANA-IV')
concat(tu3_arr,'TUCANA-III')
concat(seg_arr,'SEGUE-1')
concat(hor_arr,'HOROLOGIUM-II')
concat(boo_arr,'BOOTES-II')
