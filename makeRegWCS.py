import numpy as np
from getJdan import getJdan

workDir='./catRawMags1305/catDir/'

wcsRA, wcsDEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xc1, yc1, xt1, yt1, xr2, yr2, xc2, yc2, xo2, yo2, xr3, yr3, xc3, yc3, xo3, yo3, xr4, yr4, xc4, yc4, xo4, yo4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

wra1c, wdec1c, wra1o, wdec1o, xt2, yt2, wra2c, wdec2c, wra2t, wdec2t, xt3, yt3, wra3c, wdec3c, wra3t, wdec3t, xt4, yt4, wra4c, wdec4c, wra4t, wdec4t = 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62

wra2o, wdec2o, wra3o, wdec3o, wra4o, wdec4o = 63, 64, 65, 66, 67, 68


def regMake(targname,filt,workDir):
    jdanUse = getJdan(targname,filt)

    master = np.genfromtxt(workDir+targname+"_"+filt+"_at_long_DRCwcs.dat")

    np.savetxt(workDir+targname+'_mastCoords_'+filt+'.reg', master[:,[wcsRA,wcsDEC]],fmt='%1.6f')

    for jj, jdan in enumerate(jdanUse):

        if jj==0:
            wrc_idx = wra1c
            wdc_idx = wdec1c
            wro_idx = wra1o
            wdo_idx = wdec1o
            wrt_idx = wcsRA
            wdt_idx = wcsDEC
        elif jj==1:
            wrc_idx = wra2c
            wdc_idx = wdec2c
            wro_idx = wra2o
            wdo_idx = wdec2o
            wrt_idx = wra2t
            wdt_idx = wdec2t
        elif jj==2:
            wrc_idx = wra3c
            wdc_idx = wdec3c
            wro_idx = wra3o
            wdo_idx = wdec3o
            wrt_idx = wra3t
            wdt_idx = wdec3t
        elif jj==3:
            wrc_idx = wra4c
            wdc_idx = wdec4c
            wro_idx = wra4o
            wdo_idx = wdec4o
            wrt_idx = wra4t
            wdt_idx = wdec4t

        # dcFile = open(workDir+jdan+'_dc.reg')
        # ocFile = open(workDir+jdan+'_oc.reg')

        np.savetxt(workDir+jdan+'_dcCoords_drcWCS.reg', master[:,[wrc_idx,wdc_idx]],fmt='%1.6f')

        np.savetxt(workDir+jdan+'_ocCoords_drcWCS.reg', master[:,[wro_idx,wdo_idx]],fmt='%1.6f')

        np.savetxt(workDir+jdan+'_ltCoords_drcWCS.reg', master[:,[wrt_idx,wdt_idx]],fmt='%1.6f')

    return None

filt_arr = ['F606W', 'F814W']
targname = 'HOROLOGIUM-I'
for ff, filt in enumerate(filt_arr):
    regMake(targname,filt,workDir)
