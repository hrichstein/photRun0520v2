import numpy as np

# Constants
EEVband = 0.839 #4 pix rad
EEIband  = 0.830 #4 pix rad

ZPV = 26.667
ZPI = 26.779

seDir = '/Volumes/Spare Data/photRun0520/sDRC_2606/'

def flux2mag(flux_arr,filt):

    # Takes flux from each filter and converts to STMAG mag

    if filt=='F606W':
        EE = EEVband
        ZP = ZPV
    else:
        EE = EEIband
        ZP = ZPI

    mags = -2.5*np.log10(flux_arr/EE) + ZP


    return mags

def get_mags(targname,filt):

    if filt=='F606W':
        sex_f = seDir + "se_jdan20010_HOR-I_F606W.dat"
    else:
        sex_f = seDir + "se_jdan21010_HOR-I_F814W.dat"

    # flags,RA,DEC,x,y,flux,class_star
    cat_orig = np.loadtxt(sex_f,comments='#',dtype=None,usecols=(0,1,2,3,4,5,9,14))

    flags = np.array(cat_orig[:,0],dtype=int)
    cat = cat_orig[flags < 16] # doing a quality cut

    flux = cat[:,5]
    flux_idx = np.where(flux>=0)[0]
    flux = flux[flux_idx]

    cat = cat[flux_idx]

    newCol = np.zeros((len(cat),2)) # tacking on magnitude
                                    # and an ID
    newCol[:,0] = flux2mag(flux,filt)
    newCol[:,1] = np.arange(0,len(cat),1,dtype=int)

    cat = np.hstack((cat,newCol))

    header = "flags RA DEC xr yr flux bkgd c_star magr id"

    np.savetxt(seDir+'cat_'+targname+'_'+filt+'_'+'wMagBkgd.dat',
              cat,fmt="%d %1.7f %1.7f %1.7f %1.4f %1.4f %1.6e %1.3f %1.4f %d",header=header)


    return None

targname='HOROLOGIUM-I'
get_mags(targname,'F606W')
get_mags(targname,'F814W')




#
