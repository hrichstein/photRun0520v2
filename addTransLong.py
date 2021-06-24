from astropy.io import fits
from astropy import wcs
import numpy as np
from getJdan import getJdan

workDir = './'
upperDir = "/Volumes/Spare Data/Hannah_Data/"

filt_arr = ['F606W', 'F814W']
targname = 'HOROLOGIUM-I'

date = '1305'
catDir = workDir+'catRawMags'+date+'/catDir/'

wcsRA, wcsDEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xc1, yc1, xt1, yt1, xr2, yr2, xc2, yc2, xo2, yo2, xr3, yr3, xc3, yc3, xo3, yo3, xr4, yr4, xc4, yc4, xo4, yo4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

wra1c, wdec1c, wra1o, wdec1o, xt2, yt2, wra2c, wdec2c, wra2t, wdec2t, xt3, yt3, wra3c, wdec3c, wra3t, wdec3t, xt4, yt4, wra4c, wdec4c, wra4t, wdec4t = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21

wra2o, wdec2o, wra3o, wdec3o, wra4o, wdec4o = 22, 23, 24, 25, 26, 27

def addTranscols(targname,filt,iter,workDir='./'):

    jdanUse = getJdan(targname,filt)

    if filt=='F606W':
        fils = 'f606w/'
    elif filt=='F814W':
        fils = 'f814w/'

    cat = np.loadtxt(catDir+targname+'_'+filt+'_magList{0:d}.dat'.format(iter))

    newCol = np.zeros((len(cat),31))

    image = fits.open('/Volumes/Spare Data/photRun0520/hor1DRCs/jdan21010_drc.fits')
    w = wcs.WCS(header=image[1].header,fobj=image)

    newCol[:,0],newCol[:,1] = w.wcs_pix2world(cat[:,xc1],cat[:,yc1],1)

    # xt1,yt1 is the same as xo1 and yo1, since there's no offset for the first dither
    newCol[:,2],newCol[:,3] = w.wcs_pix2world(cat[:,xt1],cat[:,yt1],1)

    for jj in range(len(jdanUse)-1):

        if jj==0:
            xc_idx = xc2
            yc_idx = yc2
            xt_idx = xt2
            yt_idx = yt2
            xo_idx = xo2
            yo_idx = yo2
            wrc_idx = wra2c
            wdc_idx = wdec2c
            wrt_idx = wra2t
            wdt_idx = wdec2t
            wro_idx = wra2o
            wdo_idx = wdec2o
        elif jj==1:
            xc_idx = xc3
            yc_idx = yc3
            xt_idx = xt3
            yt_idx = yt3
            xo_idx = xo3
            yo_idx = yo3
            wrc_idx = wra3c
            wdc_idx = wdec3c
            wrt_idx = wra3t
            wdt_idx = wdec3t
            wro_idx = wra3o
            wdo_idx = wdec3o
        elif jj==2:
            xc_idx = xc4
            yc_idx = yc4
            xt_idx = xt4
            yt_idx = yt4
            xo_idx = xo4
            yo_idx = yo4
            wrc_idx = wra4c
            wdc_idx = wdec4c
            wrt_idx = wra4t
            wdt_idx = wdec4t
            wro_idx = wra4o
            wdo_idx = wdec4o

        transCat = np.loadtxt(workDir+catDir+jdanUse[jj+1]+"_"+targname+"_"+filt+"_t.dat")

        newCol[:,xt_idx] = transCat[:,0]
        newCol[:,yt_idx] = transCat[:,1]

        newCol[:,wrt_idx], newCol[:,wdt_idx] = w.wcs_pix2world(transCat[:,0],transCat[:,1],1)
        newCol[:,wrc_idx],newCol[:,wdc_idx] = w.wcs_pix2world(cat[:,xc_idx],cat[:,yc_idx],1)
        newCol[:,wro_idx],newCol[:,wdo_idx] = w.wcs_pix2world(cat[:,xo_idx],cat[:,yo_idx],1)

        image.close()

    print(newCol[0][wra2o])

    newCol[:,-3] = np.mean( np.array([newCol[:,wra2t],newCol[:,wra3t],newCol[:,wra4t]]))
    newCol[:,-2] = np.mean( np.array([newCol[:,wdec2t],newCol[:,wdec3t],newCol[:,wdec4t]]))
    newCol[:,-1] = np.nanmedian(np.array([cat[:,mag1],cat[:,mag2],cat[:,mag3],cat[:,mag4] ]))

    print(newCol[0])

    cat = np.hstack((cat, newCol))


    header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1 xr2 yr2 xc2 yc2 xo2 yo2 xr3 yr3 xc3 yc3 xo3 yo3 xr4 yr4 xc4 yc4 xo4 yo4 wra1c wdec1c wra1o wdec1o xt2 yt2 wra2c wdec2c wra2t wdec2t xt3 yt3 wra3c wdec3c wra3t wdec3t xt4 yt4 wra4c wdec4c wra4t wdec4t wra2o wdec2o wra3o wdec3o wra4o wdec4o meanRA_234 meanDEC_234 median_mag'


    # print(cat[0][wra2c+41])

    np.savetxt(workDir+catDir+targname+"_"+filt+"_at_long_DRCwcs.dat", cat,header=header)


    return None

for ff, filt in enumerate(filt_arr):
    addTranscols(targname,filt,iter=1,workDir='./')
