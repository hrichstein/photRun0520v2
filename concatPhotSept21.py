import numpy as np
import matplotlib.pyplot as plt

# def f2mag_dirs(targname,date='20Aug',workDir='./'):
#
#     return workDir+'catRawMags'+date+'/catDir_'+targname+'/'
catDir = 'catsSept21/'

outDir = './catsSept21/'

tri_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']
tuc_arr = ['TUCANA-II-NE','TUCANA-II-SE']

def concat(name_arr,source):

    for ff in range(len(name_arr)):
        # seDir, magCatDir, catDir = f2mag_dirs(name_arr[ff],date='20Aug',workDir='./')
        if ff==0:
            f_out = np.genfromtxt(catDir+'catsSept21_'+name_arr[ff]+'.dat')
            fhead = np.genfromtxt(catDir+'catsSept21_'+name_arr[ff]+'.dat',names=True)
        else:
            f_temp = np.genfromtxt(catDir+'catsSept21_'+name_arr[ff]+'.dat')
            f_out = np.vstack((f_out,f_temp))

    colNs = np.array(fhead.dtype.names)
    header = ' '.join(colNs)


    np.savetxt(outDir+source+"_all_source.dat",f_out,header=header,fmt='%1.5f')

    makeCMD(outDir,source)

    return None


def makeCMD(dir,source):

    file = np.genfromtxt(dir+source+"_all_source.dat",names=True)

    f_idx = np.logical_and(file['c_star_f606w']>=0.5, file['c_star_f814w']>=0.5)

    flc_g = file[f_idx]

    fig, ax = plt.subplots(figsize=(4,6.5))

    ax.scatter(flc_g['magZPT_f606w']-flc_g['magZPT_f814w'],flc_g['magZPT_f606w'],s=5,label='FLC')

    ax.set_ylim(28,17.5)
    ax.set_xlim(-1.5,1.5)
    ax.set_xlabel('F606W-F814W')
    ax.set_ylabel('F606W')

    ax.set_title(source)

    ax.legend()

    plt.savefig(dir+source+'_all_cmd.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None

concat(tri_arr,'TRIANGULUM-II')
concat(tuc_arr,'TUCANA-II')
