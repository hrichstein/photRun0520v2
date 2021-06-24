import numpy as np

def pRefStars(targname,filt,iter,catDir='./',magHi=24.5,magLo=21,\
    stdTol=1,posTol=1):

    data = np.genfromtxt(catDir+targname+'-cut-stdPixMatched_'+filt+'_{0:d}.dat'.format(iter-1),names=True)

    keep = data['mean']!=data['mean']
    for ll in range(len(data)):
        temp_keep1 = np.logical_and(np.logical_and(data['mean']>=magLo,\
                data['mean']<=magHi),np.logical_and(data['pix_std']<posTol,data['stdev']<stdTol))

        keep = np.logical_or(keep,temp_keep1)

    print('Number of Reference Stars:',len(data[keep]))

    header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1  xr2 yr2 xc2 yc2 xt2 yt2 xr3 yr3 xc3 yc3 xt3 yt3 xr4 yr4 xc4 xc4 xt4 yt4 mean stdev pix_std coord_std cut_flag idx_cut num_abv_std'

    form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %d %d %d'

    np.savetxt(catDir+targname+'_refStars_'+filt+'_{0:d}.dat'.format(iter), data[keep], fmt=form, header=header)

    return None
