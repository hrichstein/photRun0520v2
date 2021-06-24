import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import vstack
from photutils import aperture_photometry,CircularAnnulus,CircularAperture
from photutils import DAOStarFinder
from drizzlepac import pixtopix

# flcDir will be ../ for where this code is


def runPhotUtils(targname,filt,radius=4,jdanUse=None,drcFile='None.dat',
                 saveDir='./',flcDir='./',suffix='_pu.dat'):

    if filt=='F606W':
        if abs(radius - 4) <= 1e-3:
            EEband = 0.839  # 4 pixel
        elif abs(radius - 3) <= 1e-3:
            EEband = 0.795
        ZPT = 26.667
        fils = '_f606w/'

    elif filt=='F814W':
        if abs(radius - 4) <= 1e-3:
            EEband = 0.830  # 4 pixels
        elif abs(radius - 3) <= 1e-3:
            EEband = 0.77
        ZPT = 26.779
        fils = '_f814w/'

    for ff in range(len(jdanUse)):
        image = flcDir + targname + fils + '/crClean/' + jdanUse[ff] \
            + '_crclean.fits'

        hdu = fits.open(image)
        sci1 = hdu[1].data  # chip 2
        sci2 = hdu[4].data  # chip 1
        hdr = hdu[0].header
        hdu.close()

        sci_arr = [sci1,sci2]

        for cc, sci in enumerate(sci_arr):
            data = sci.copy()

            mean, median, std = sigma_clipped_stats(data, sigma=3.0,
                                                    maxiters=10)

            daofind = DAOStarFinder(fwhm=2., threshold=5. * std)
            sources = daofind(data - median)

            loc = np.array([sources['xcentroid'], sources['ycentroid']])
            positions = np.transpose(loc)

            apertures_rad = CircularAperture(positions, r=4.)
            rawflux_rad = aperture_photometry(data, apertures_rad)

            annulus_apertures = CircularAnnulus(positions, r_in=9.,
                                                r_out=12.)

            annulus_masks = annulus_apertures.to_mask(method='center')

            bkg_median = []
            for mask in annulus_masks:

                annulus_data = mask.multiply(data)
                annulus_data_1d = annulus_data[mask.data > 0]
                _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
                bkg_median.append(median_sigclip)

            bkg_median = np.array(bkg_median)

            rawflux_rad['annulus_median'] = bkg_median
            rawflux_rad['aper_bkg'] = bkg_median * apertures_rad.area

            rawflux_rad['final_phot'] = rawflux_rad['aperture_sum'] \
                - rawflux_rad['aper_bkg']

            mask_negative = (rawflux_rad['final_phot'] > 0)
            rawflux_pos_rad = rawflux_rad[mask_negative]

            final_phot = -2.5 * np.log10(rawflux_pos_rad['final_phot']
                                         / EEband) \
                + 2.5 * np.log10(hdr['exptime']) + ZPT

            xin = np.zeros((len(rawflux_pos_rad),1))
            xin[:,0] = rawflux_pos_rad['xcenter']

            yin = np.zeros((len(rawflux_pos_rad),1))
            yin[:,0] = rawflux_pos_rad['ycenter']
            xt1,yt1 = pixtopix.tran(flcDir + targname + fils + '/crClean/'
                                    + jdanUse[ff]
                                    + '_crclean.fits[sci,{0:d}]'.format(cc+1),
                                    drcFile + '[sci,1]','forward',
                                    x=xin.flatten(), y=yin.flatten())

            rawflux_pos_rad['xDRC'] = xt1
            rawflux_pos_rad['yDRC'] = yt1

            rawflux_pos_rad['magr'] = final_phot

            if cc==0:
                temp_store = rawflux_pos_rad
            else:
                out_arr = vstack([temp_store,rawflux_pos_rad])

        out_arr['id'] = np.arange(0,len(out_arr),1)

        s0 = ' '
        header = s0.join(out_arr.dtype.names)

        outname = saveDir + jdanUse[ff] + '_' + filt + suffix

        np.savetxt(outname,out_arr,header=header)

    return None


# def f2mag_dirs(targname,date='10Aug',workDir='./'):
#
#     # Gives the name of the directory with se files
#     # Creates folder for raw mags if it doesn't exists
#     # returns the names of both
#
#     magCatDir = workDir + 'photUtils' + date + '/'
#     catDir = magCatDir + 'catDir_' + targname + '/'
#
#     if not os.path.exists(os.path.join(".",magCatDir)):
#         os.makedirs(magCatDir)
#     if not os.path.exists(os.path.join(".",catDir)):
#         os.makedirs(catDir)
#
#
#     return catDir

# targname = 'HOROLOGIUM-I'
# jdanUse_814 = getJdan(targname,filt='F814W')
# jdanUse_606 = getJdan(targname,filt='F606W')
#
# for jj in range(len(jdanUse_814)):
#     jdan = jdanUse_814[jj]
#     runPhotUtils(targname,jdan,'F814W')
#
# for jj in range(len(jdanUse_606)):
#     jdan = jdanUse_606[jj]
#     runPhotUtils(targname,jdan,'F606W')
