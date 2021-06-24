import numpy as np
import os
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry,CircularAnnulus,CircularAperture
from photutils import DAOStarFinder


# Work in progress
# Assuming 4 pixel aperture

drcDir = '../Hannah_Data/origDRCs/'

def runPhotUtils(drcInfo,saveDir='./'):

    info = np.loadtxt(drcInfo,dtype=str)
    infoN = np.genfromtxt(drcInfo,dtype=str,names=True)
    nameCols = np.array(infoN.dtype.names)

    jdan = np.int(np.where(nameCols=='JDAN')[0])
    filt1 = np.int(np.where(nameCols=='FILTER1')[0])
    filt2 = np.int(np.where(nameCols=='FILTER2')[0])
    targN = np.int(np.where(nameCols=='TARGNAME')[0])

    fileNames = info[:,jdan]

    for ff in range(len(info)):

        image = drcDir + fileNames[ff] + ".fits"
        print(image)

        hdu = fits.open(image)
        sci = hdu[1].data
        hdr = hdu[1].header
        hdu.close()

        data = sci.copy()

        f1 = info[:,filt1][ff]
        f2 = info[:,filt2][ff]

        targname = info[:,targN][ff]
        print(targname)

        saveDir = f2mag_dirs(targname,date='20Aug',workDir='./')

        if (f1=='F606W') or (f2=='F606W'):
            filt = 'F606W'
            EEband = 0.839
            ZPT = 26.667

        elif (f1=='F814W') or (f2=='F814W'):
            filt = 'F814W'
            EEband = 0.830
            ZPT = 26.779

        mean, median, std = sigma_clipped_stats(data, sigma=3.0, \
                                                maxiters=10)

        daofind = DAOStarFinder(fwhm=2.5, threshold=5.*std)
        sources = daofind(data - median)

        loc = np.array([sources['xcentroid'], sources['ycentroid']])
        positions = np.transpose(loc)

        apertures_r4 = CircularAperture(positions, r=4.)
        rawflux_r4 = aperture_photometry(data, apertures_r4)

        annulus_apertures = CircularAnnulus(positions, r_in=9., \
                                            r_out=12.)

        annulus_masks = annulus_apertures.to_mask(method='center')

        bkg_median = []
        for mask in annulus_masks:

            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)

        bkg_median = np.array(bkg_median)

        rawflux_r4['annulus_median'] = bkg_median
        rawflux_r4['aper_bkg'] = bkg_median*apertures_r4.area

        rawflux_r4['final_phot'] = rawflux_r4['aperture_sum'] \
        - rawflux_r4['aper_bkg']

        mask_negative = (rawflux_r4['final_phot'] > 0)
        rawflux_pos_r4 = rawflux_r4[mask_negative]

        final_phot = -2.5*np.log10(rawflux_pos_r4['final_phot']/EEband) + ZPT

        rawflux_pos_r4['magr'] = final_phot

        rawflux_pos_r4['id'] = np.arange(0,len(rawflux_pos_r4),1)

        s0 = ' '
        header = s0.join(rawflux_pos_r4.dtype.names)

        outname = saveDir + fileNames[ff] + '_' + filt + 'photU.dat'

        np.savetxt(outname,rawflux_pos_r4,header=header)

        print('Moving On.')


    return None


def f2mag_dirs(targname,date='20Aug',workDir='./'):

    # Gives the name of the directory with se files
    # Creates folder for raw mags if it doesn't exists
    # returns the names of both

    magCatDir = workDir + 'photUtils' + date + '/'
    catDir = magCatDir + 'catDir_' + targname + '/'

    if not os.path.exists(os.path.join(".",magCatDir)):
        os.makedirs(magCatDir)
    if not os.path.exists(os.path.join(".",catDir)):
        os.makedirs(catDir)


    return catDir


# targname_arr = np.genfromtxt('targnamesDirections.txt',dtype='str')

targname_arr = np.genfromtxt('targnamesPost.txt',dtype='str')
#
def main():
#
    for c1,targname in enumerate(targname_arr):
        saveDir = f2mag_dirs(targname,date='28Sep',workDir='./')
#
    runPhotUtils("../Hannah_Data/drcTargInfo2.dat",saveDir='./')


if __name__ == '__main__':
    main()
