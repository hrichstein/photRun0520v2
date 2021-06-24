import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry,CircularAnnulus,CircularAperture
from photutils import DAOStarFinder

from getJdan import getJdan

def runPhotUtils(targname,jdan,filt):

    if filt=='F606W':
        EEband = 0.839
        ZPT = 26.667
        fils = '_f606w/'

    elif filt=='F814W':
        EEband = 0.830
        ZPT = 26.779
        fils = '_f814w/'

    image = targname + fils + '/crClean/' + jdan + '_WJ2.fits'

    hdu = fits.open(image)
    sci = hdu[0].data
    hdr = hdu[0].header
    hdu.close()

    data = sci.copy()

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

    final_phot = -2.5*np.log10(rawflux_pos_r4['final_phot']/EEband) \
        +2.5*np.log10(hdr['exptime']) + ZPT

    rawflux_pos_r4['magr'] = final_phot
    s0 = ' '
    header = s0.join(rawflux_pos_r4.dtype.names)

    saveDir = 'photUtils0820/hor1run0608/'
    outname = saveDir + jdan + '_' + filt + 'photU.dat'

    np.savetxt(outname,rawflux_pos_r4,header=header)


    return None


targname = 'HOROLOGIUM-I'
jdanUse_814 = getJdan(targname,filt='F814W')
jdanUse_606 = getJdan(targname,filt='F606W')

for jj in range(len(jdanUse_814)):
    jdan = jdanUse_814[jj]
    runPhotUtils(targname,jdan,'F814W')

for jj in range(len(jdanUse_606)):
    jdan = jdanUse_606[jj]
    runPhotUtils(targname,jdan,'F606W')
