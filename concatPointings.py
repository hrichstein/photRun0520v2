import numpy as np
import matplotlib.pyplot as plt
from f2mag0707 import f2mag_dirs

outDir = './catDir0720/'

tri_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']
tuc_arr = ['TUCANA-II-NE',
'TUCANA-II-NW','TUCANA-II-SE','TUCANA-II-SW']

def concat(name_arr,source):

    for ff in range(len(name_arr)):
        seDir, magCatDir, catDir = f2mag_dirs(name_arr[ff],date='1305',workDir='./')
        if ff==0:
            f_out = np.genfromtxt(catDir+name_arr[ff]+'_sourceList0720.dat')
        else:
            f_temp = np.genfromtxt(catDir+name_arr[ff]+'_sourceList0720.dat')
            f_out = np.vstack((f_out,f_temp))

    header = 'RA_f606w DEC_f606w flags_f606w c_star_f606w xt1_f606w yt1_f606w mean_f606w stdev_f606w magZPT_f606w magZPTerr_f606w RA_f814w DEC_f814w flags_f814w c_star_f814w xt1_f814w yt1_f814w mean_f814w stdev_f814w magZPT_f814w magZPTerr_f814w'

    form = "%1.7f %1.7f %d %1.3f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.7f %1.7f %d %1.3f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f"

    np.savetxt(outDir+source+"_all_sourceList0720.dat",f_out,header=header,fmt=form)

    makeCMD(outDir,source)

    return None


def makeCMD(dir,source):

    file = np.genfromtxt(dir+source+"_all_sourceList0720.dat",names=True)

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
