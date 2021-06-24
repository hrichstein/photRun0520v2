from astropy.io import fits
from astropy import wcs
import numpy as np
import os

from hst_func import *

def distCor(targname,filt,workDir='./'):

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):
        #Load catalog
        cat = np.genfromtxt(workDir+ jdanUse[ff]+'_'+targname+'_'+filt+'_wMag.dat',names=True)
        catCat = np.loadtxt(workDir+ jdanUse[ff]+'_'+targname+'_'+filt+'_wMag.dat')
        #Do correction
        cor = acsDistortion(upperDir+'wfc_'+filt,cat['xr'],cat['yr'])
        #Add columns
        cat = np.hstack((catCat,cor))

        header = "flags RA DEC xr yr flux c_star magr id xc yc"
        form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f"

        np.savetxt(workDir+jdanUse[ff]+'_'+targname+'_'+filt+'_dc.dat',cat,fmt=form,header=header)

    return None

def offCor(targname,filt,workDir='./'):

    jdanUse = getJdan(targname,filt)

    if filt=='F606W':
        fils = 'f606w/'
    elif filt=='F814W':
        fils = 'f814w/'

    for jj, jdan in enumerate(jdanUse):

        # Need to reference whatever directory you have the flc fits files in
            # or copy the flcs to this directory
        # Load images, retrieve offset info
        tempim = fits.open(workDir+targname+'_'+fils+jdan+"_flc.fits")

        xoff = float(tempim[0].header["POSTARG1"])
        yoff = float(tempim[0].header["POSTARG2"])

        # Load the respective catalog
        cat = np.genfromtxt(workDir+jdan+"_"+targname+"_"+filt+"_dc.dat",names=True)
        catCat = np.loadtxt(workDir+jdan+"_"+targname+"_"+filt+"_dc.dat")

        # Create an array for the new values
        newCol = np.zeros((len(cat),2))

        # Apply offsets to columns
        newCol[:,0] = cat['xc'] - (offset * xoff)
        newCol[:,1] = cat['yc'] - (offset * yoff)

        # Combine to single array and save out
        cat = np.hstack((catCat, newCol))

        header = "flags RA DEC xr yr flux c_star magr id xc yc xo yo"
        form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f"

        np.savetxt(workDir+jdan+"_"+targname+"_"+filt+"_oc.dat",cat,header=header,fmt=form)


    return None


distCor('HOROLOGIUM-I','F814W',workDir='./')
offCor('HOROLOGIUM-I','F814W',workDir='./')
