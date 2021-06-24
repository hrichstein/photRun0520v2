from astropy.io import fits
from astropy import wcs
import numpy as np

workDir='./catRawMags1305/catDir/'
drcDir = '/Volumes/Spare Data/photRun0520/hor1DRCs/'
file = 'jdan21010_drc.fits'

dir3  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
drc_file = np.genfromtxt(dir3 + 'HOROLOGIUM-I_sfErr.dat')

RA_v, DEC_v, x_v, y_v, fAper_v, fErr_v, magAper_v, magErr_v, magRaw_v, magRed_v, magAbs_v, elong_v, ellip_v, class_Star_v, RA_i, DEC_i, x_i, y_i, fAper_i, fErr_i, magAper_i, magErr_i, magRaw_i, magRed_i, magAbs_i, elong_i, ellip_i, class_Star_i, corrF_errV, corrF_errI, corrM_errV, corrM_errI = 0, 1, 2 ,3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31


image = fits.open(drcDir+file)
w = wcs.WCS(header=image[1].header,fobj=image)


newCols = np.zeros((len(drc_file),2))

newCols[:,0], newCols[:,1] = w.wcs_pix2world(drc_file[:,x_i],drc_file[:,y_i],1)

image.close()

drcidx = np.argsort(drc_file[:,magRaw_i])[:1000]
drc_1000 = newCols[drcidx]

np.savetxt(workDir+'HOROLOGIUM-I_F814W_seDRCwcs.reg',drc_1000,fmt='%1.6f')
