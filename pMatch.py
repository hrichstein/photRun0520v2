import numpy as np

try:
    from scipy.spatial import cKDTree as KDT
except ImportError:
    from scipy.spatial import KDTree as KDT

upperDir = "/Volumes/Spare Data/Hannah_Data/"
date = '1305'

wcsRA, wcsDEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xc1, yc1, xt1, yt1, xr2, yr2, xc2, yc2, xt2, yt2, xr3, yr3, xc3, yc3, xt3, yt3, xr4, yr4, xc4, yc4, xt4, yt4, mean, stdev, pix_std, coord_std, cut_flag, idx_cut, num_abv_std  = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47

def pix_match(targname,iter,catDir='./',matchtol=2):

    infileF606W = catDir+targname+'_F606W_cut_std_{0:d}.dat'.format(iter)
    infileF814W = catDir+targname+'_F814W_cut_std_{0:d}.dat'.format(iter)

    f606w = np.loadtxt(infileF606W, dtype='float')
    f814w = np.loadtxt(infileF814W, dtype='float')

    print('Length F606W',len(f606w))
    print('Length F814W',len(f814w))

    xt_f606w = np.vstack((f606w[:,xt1],f606w[:,xt2],f606w[:,xt3],f606w[:,xt4]))
    xt_f814w = np.vstack((f814w[:,xt1],f814w[:,xt2],f814w[:,xt3],f814w[:,xt4]))

    yt_f606w = np.vstack((f606w[:,yt1],f606w[:,yt2],f606w[:,yt3],f606w[:,yt4]))
    yt_f814w = np.vstack((f814w[:,yt1],f814w[:,yt2],f814w[:,yt3],f814w[:,yt4]))

    xm_606 = np.nanmean(xt_f606w,axis=0)
    xm_814 = np.nanmean(xt_f814w,axis=0)

    ym_606 = np.nanmean(yt_f606w,axis=0)
    ym_814 = np.nanmean(yt_f814w,axis=0)

    coords1 = np.empty((xm_606.size,2))
    coords2 = np.empty((xm_814.size,2))

    coords1[:,0] = xm_606
    coords1[:,1] = ym_606

    coords2[:,0] = xm_814
    coords2[:,1] = ym_814

    kdt = KDT(coords1)
    idxs2 = kdt.query(coords2)[1]

    ds = distArr(xm_814,ym_814,xm_606[idxs2],ym_606[idxs2])

    idxs1 = np.arange(xm_814.size)

    msk = ds < matchtol
    idxs1 = idxs1[msk]
    idxs2 = idxs2[msk]
    ds = ds[msk]

    # print(max(idxs1))
    # print(max(idxs2))

    outfile = catDir+targname+'-cut_F606W_match_{0:d}.dat'.format(iter)
    np.savetxt(outfile, idxs2, fmt='%4i')

    outfile = catDir+targname+'-cut_F814W_match_{0:d}.dat'.format(iter)
    np.savetxt(outfile, idxs1, fmt='%4i')

    outfile = catDir+targname+'-cut_ds_{0:d}.dat'.format(iter)
    np.savetxt(outfile, ds, fmt='%1.4f')


    return None

def matchMags(targname,iter,catDir='./'):

    filters = ['F606W','F814W']

    header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1 xr2 yr2 xc2 yc2 xt2 yt2 xr3 yr3 xc3 yc3 xt3 yt3 xr4 yr4 xc4 xc4 xt4 yt4 mean stdev pix_std coord_std cut_flag idx_cut num_abv_std'

    form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %d %d %d'

    for ff in range(len(filters)):
        dataFile = np.genfromtxt(catDir+targname+"_"+filters[ff]+'_cut_std_{0:d}.dat'.format(iter))
        idxFile = np.genfromtxt(catDir+targname+'-cut_'+filters[ff]+"_match_{0:d}.dat".format(iter),dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(catDir+targname+'-cut-stdPixMatched_'+filters[ff]+'_{0:d}.dat'.format(iter),outArr,fmt=form,header=header)

    return None


def distArr(x0,y0,x_arr,y_arr):
    dist_arr = np.sqrt( (x0-x_arr)**2 + (y0-y_arr)**2 )

    return dist_arr


def pMatch(targname,iter,catDir='./',matchtol=2):

    pix_match(targname,iter,catDir,matchtol)
    matchMags(targname,iter,catDir)

    return None
