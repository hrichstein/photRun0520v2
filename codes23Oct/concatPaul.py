import numpy as np
import matplotlib.pyplot as plt

# def f2mag_dirs(targname,date='20Aug',workDir='./'):
#
#     return workDir+'catRawMags'+date+'/catDir_'+targname+'/'
catDir = 'drcPhot02Nov/'

# plotDir = catDir + 'paulCMDs/'

tri_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']
tuc_arr = ['TUCANA-II-NE','TUCANA-II-SE']
urs_arr = ['URSA-MAJOR-II-WEST','URSA-MAJOR-II-EAST']
tu4_arr = ['TUCANA-IV-SOUTH','TUCANA-IV-NORTH']
tu3_arr = ['TUCANA-III-WEST','TUCANA-III-EAST']
seg_arr = ['SEGUE-1-WEST','SEGUE-1-EAST']
hor_arr = ['HOROLOGIUM-II-WEST','HOROLOGIUM-II-EAST']
boo_arr = ['BOOTES-II-SOUTH','BOOTES-II-NORTH']


def concat(name_arr,source,catDir='../drcPhot02Nov/'):

    for ff in range(len(name_arr)):

        if ff==0:
            f_out = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                  + name_arr[ff] + '_matchedDRCfilt.dat')
            fhead = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                  + name_arr[ff] + '_matchedDRCfilt.dat',
                                  names=True)
        else:
            f_temp = np.genfromtxt(catDir+ 'catDir_' + name_arr[ff] + '/'
                                   + name_arr[ff] + '_matchedDRCfilt.dat')
            f_out = np.vstack((f_out,f_temp))

    colNs = np.array(fhead.dtype.names)
    header = ' '.join(colNs)

    np.savetxt(catDir+ 'catDir_' + name_arr[ff] + '/' + source
               + "_matchedDRCfilt_ALL.dat",f_out,header=header,fmt='%1.5f')

    makeCMD(catDir+ 'catDir_' + name_arr[ff] + '/',source)

    return None


def makeCMD(dir,source):

    file = np.genfromtxt(dir + source + "_matchedDRCfilt_ALL.dat",names=True)

    fig, ax = plt.subplots(figsize=(4,6.5))

    ax.scatter(file['magr_f606w']-file['magr_f814w'],file['magr_f606w'],s=5)

    ax.set_ylim(28,17.5)
    ax.set_xlim(-1.5,1.5)
    ax.set_xlabel('F606W-F814W')
    ax.set_ylabel('F606W')

    ax.set_title(source)

    plt.savefig('../drcPhot02Nov/paulCMDs/'+source+'_drc_ALL.png',dpi=600,
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
