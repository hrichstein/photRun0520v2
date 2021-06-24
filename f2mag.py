import numpy as np
from getJdan import getJdan
import os

# workDir = './'
# upperDir = "/Volumes/Spare Data/Hannah_Data/"

#   1 FLAGS                  Extraction flags
#   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#   4 X_IMAGE                Object position along x                                    [pixel]
#   5 Y_IMAGE                Object position along y                                    [pixel]
#   6 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
#   7 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
#   8 MAG_APER               Fixed aperture magnitude vector                            [mag]
#   9 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#  10 BACKGROUND             Background at centroid position                            [count]
#  11 THRESHOLD              Detection threshold above background                       [count]
#  12 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
#  13 ELONGATION             A_IMAGE/B_IMAGE
#  14 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
#  15 CLASS_STAR             S/G classifier output

# Constants
EEVband = 0.795 #3 pix rad
EEIband  = 0.770 #3 pix rad

ZPV = 26.667
ZPI = 26.779

def get_mags(targname,filt,date,workDir='./'):

    jdanUse = getJdan(targname,filt)
    seDir, saveDir, catDir = f2mag_dirs(date=date,workDir=workDir)

    for ff in range(len(jdanUse)):
        #Load sextractor file catalog
        sex_f = workDir+seDir+"se_"+jdanUse[ff]+'_'+targname+'_'+filt+'.dat'

        # flags,RA,DEC,x,y,flux,class_star
        cat_orig = np.loadtxt(sex_f,comments='#',dtype=None,usecols=(0,1,2,3,4,5,14))

        flags = np.array(cat_orig[:,0],dtype=int)
        cat = cat_orig[flags < 16] # doing a quality cut

        flux = cat[:,5]
        exptime = get_expt(jdanUse[ff],filt,workDir=workDir)

        newCol = np.zeros((len(cat),2)) # tacking on magnitude
                                        # and an ID
        newCol[:,0] = flux2mag(flux,filt,exptime)
        newCol[:,1] = np.arange(0,len(cat),1,dtype=int)

        cat = np.hstack((cat,newCol))

        header = "flags RA DEC xr yr flux c_star magr id"

        np.savetxt(workDir+saveDir+jdanUse[ff]+'_'+targname+'_'+filt+'_'+'wMag.dat',
                  cat,fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d",header=header)


    return None


def f2mag_dirs(date='1305',workDir='./'):

    # Gives the name of the directory with se files
    # Creates folder for raw mags if it doesn't exists
    # returns the names of both

    seDir = workDir + 'sFLCf' + date + '/'
    magCatDir = workDir + 'catRawMags' + date + '/'
    catDir = magCatDir + 'catDir/'

    if not os.path.exists(os.path.join(".",magCatDir)):
        os.makedirs(magCatDir)
    if not os.path.exists(os.path.join(".",catDir)):
        os.makedirs(catDir)


    return seDir, magCatDir, catDir


def flux2mag(flux_arr,filt,exptime):

    # Takes flux from each filter and converts to STMAG mag

    if filt=='F606W':
        EE = EEVband
        ZP = ZPV
    else:
        EE = EEIband
        ZP = ZPI

    mags = -2.5*np.log10(flux_arr/EE) + 2.5*np.log10(exptime) + ZP


    return mags


def get_expt(jdan,filt,workDir='./'):

    if filt=='F606W':
        expT_file = 'WJC_F606WList.dat'
    elif filt=='F814W':
        expT_file = 'WJC_F814WList.dat'

    # flcName and EXPTIME columns
    file = np.loadtxt(workDir+expT_file,usecols=(0,4),
                         comments='#',dtype=str)
    targname = jdan+'_WJ2.fits'

    exptime = file[:,1][file[:,0]==targname]
    # gets the EXPTIME (column value) for the matching flcName (row)

    return float(exptime)
