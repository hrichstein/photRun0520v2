import numpy as np
import matplotlib.pyplot as plt

def outputElen(targname, dir='./'):

    # fileN = np.genfromtxt(dir+targname+'_allMatchedZPTed.dat',names=True)
    # file = np.genfromtxt(dir+targname+'_allMatchedZPTed.dat')
    fileN = np.genfromtxt(dir+targname+'_allMatchedZPTed_mDc.dat',names=True)
    file = np.genfromtxt(dir+targname+'_allMatchedZPTed_mDc.dat')

    colNs = np.array(fileN.dtype.names)

    ra_f606w = np.int(np.where(colNs=='RA_f606w')[0])
    dec_f606w = np.int(np.where(colNs=='DEC_f606w')[0])
    flags_f606w = np.int(np.where(colNs=='flags_f606w')[0])
    c_star_f606w = np.int(np.where(colNs=='c_star_f606w')[0])
    xt1_f606w = np.int(np.where(colNs=='xt1_f606w')[0])
    yt1_f606w = np.int(np.where(colNs=='yt1_f606w')[0])
    mean_f606w = np.int(np.where(colNs=='mean_f606w')[0])
    stdev_f606w = np.int(np.where(colNs=='stdev_f606w')[0])
    magZPT_f606w = np.int(np.where(colNs=='magZPT_f606w')[0])
    magZPTerr_f606w = np.int(np.where(colNs=='magZPTerr_f606w')[0])

    ra_f814w = np.int(np.where(colNs=='RA_f814w')[0])
    dec_f814w = np.int(np.where(colNs=='DEC_f814w')[0])
    flags_f814w = np.int(np.where(colNs=='flags_f814w')[0])
    c_star_f814w = np.int(np.where(colNs=='c_star_f814w')[0])
    xt1_f814w = np.int(np.where(colNs=='xt1_f814w')[0])
    yt1_f814w = np.int(np.where(colNs=='yt1_f814w')[0])
    mean_f814w = np.int(np.where(colNs=='mean_f814w')[0])
    stdev_f814w = np.int(np.where(colNs=='stdev_f814w')[0])
    magZPT_f814w = np.int(np.where(colNs=='magZPT_f814w')[0])
    magZPTerr_f814w = np.int(np.where(colNs=='magZPTerr_f814w')[0])

    outArr = file[:,[ra_f606w,dec_f606w,flags_f606w,c_star_f606w,xt1_f606w,yt1_f606w,mean_f606w,stdev_f606w,magZPT_f606w,magZPTerr_f606w,ra_f814w,dec_f814w,flags_f814w,c_star_f814w,xt1_f814w,yt1_f814w,mean_f814w,stdev_f814w,magZPT_f814w,magZPTerr_f814w]]

    header = 'RA_f606w DEC_f606w flags_f606w c_star_f606w xt1_f606w yt1_f606w mean_f606w stdev_f606w magZPT_f606w magZPTerr_f606w RA_f814w DEC_f814w flags_f814w c_star_f814w xt1_f814w yt1_f814w mean_f814w stdev_f814w magZPT_f814w magZPTerr_f814w'

    form = "%1.7f %1.7f %d %1.3f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.7f %1.7f %d %1.3f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f"

    # np.savetxt(dir+targname+'_sourceList0720.dat',outArr,header=header,fmt=form)
    np.savetxt(dir+targname+'_sourceList0723_mDc.dat',outArr,header=header,fmt=form)


    return None


#
