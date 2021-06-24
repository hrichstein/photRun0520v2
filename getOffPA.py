from astropy.io import fits
import numpy as np
from getJdan import getJdan



scale = 20.0 #20 pixels per arcsecond

def getOffPA(targname,filt,workDir='./'):

    jdanUse = getJdan(targname,filt)

    if filt=='F606W':
        fils = 'f606w/'
    elif filt=='F814W':
        fils = 'f814w/'

    file = open(workDir+'offsetPA_'+targname+'_'+filt+'.dat','w')
    file.write('xoff yoff pa xpix ypix\n')

    for jj, jdan in enumerate(jdanUse):

        # Load images, retrieve offset info
        tempim = fits.open(workDir+targname+'_'+fils+jdan+"_flc.fits")

        xoff = float(tempim[0].header["POSTARG1"])
        yoff = float(tempim[0].header["POSTARG2"])

        pa = float(tempim[0].header["PA_V3"])

        outline = np.array([xoff, yoff, pa, xoff*scale, yoff*scale])

        file.write('{0:f} {1:f} {2:f} {3:f} {4:f} \n'.format(outline[0], outline[1], outline[2], outline[3], outline[4]))

    file.close()

    return None

filt_arr = ['F606W','F814W']
for ff, filt in enumerate(filt_arr):
    getOffPA('HOROLOGIUM-I',filt,workDir='./')
