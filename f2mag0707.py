import numpy as np
from getJdan import getJdan
import os

EEVband = 0.795 #3 pix rad
EEIband  = 0.770 #3 pix rad

ZPV = 26.667
ZPI = 26.779

def f2mag_dirs(targname,date='1305',workDir='./'):

    # Gives the name of the directory with se files
    # Creates folder for raw mags if it doesn't exists
    # returns the names of both

    seDir = workDir + 'sFLCf' + date + '/'
    magCatDir = workDir + 'catRawMags' + date + '/'
    catDir = magCatDir + 'catDir_' + targname + '/'

    if not os.path.exists(os.path.join(".",magCatDir)):
        os.makedirs(magCatDir)
    if not os.path.exists(os.path.join(".",catDir)):
        os.makedirs(catDir)


    return seDir, magCatDir, catDir


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


def get_mags(targname,filt,date,workDir='./'):

    jdanUse = getJdan(targname,filt)
    seDir, saveDir, catDir = f2mag_dirs(targname,date=date,workDir=workDir)

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

        np.savetxt(workDir+catDir+jdanUse[ff]+'_'+targname+'_'+filt+'_'+'wMag.dat',
                  cat,fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d",header=header)


    return None
