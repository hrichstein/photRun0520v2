#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import pylab as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Ellipse
from matplotlib.mlab import rec2csv

import emcee
import corner

import astropy.wcs
from astropy.io import fits

# WCS (since ugali model works in sky coordinates)
header = dict(
NAXIS   = 2,
NAXIS1 = 7902,
NAXIS2 = 7845,
CD1_1   = -9.7207414806871E-06, 
CD1_2   = 1.69676173695811E-07, 
CD2_1   = 1.69676173695811E-07, 
CD2_2   = 9.72074148068713E-06, 
CRVAL1  =    56.09056245173978, 
CRVAL2  =   -43.53022649816447, 
CRPIX1  =    3951.096595646065, 
CRPIX2  =    3922.701986935564, 
CTYPE1  = 'RA---TAN'          , 
CTYPE2  = 'DEC--TAN'          , 
ORIENTAT=   0.9999999999999952, 
)
WCS = astropy.wcs.WCS(header)

# Gaps in the image
IMAGE = Path([ [406,2390],
              [2190,7814],
              [7662,5536],
              [5704,5],
              [406,2390]])
GAP = Path([ [4842, 6714],
             [4905, 6686],
             [3037, 1208],
             [2977, 1233],
             [4842, 6714]])
IMG_BIT = 0b01
GAP_BIT = 0b10

# Bin edges
NBINS = 267
XMIN,XMAX = 0,8010
YMIN,YMAX = 0,8010
XEDGE = np.linspace(XMIN,XMAX,NBINS)
YEDGE = np.linspace(XMIN,XMAX,NBINS)

# Bin centers
XCENT = (XEDGE[1:] + XEDGE[:-1])/2.
YCENT = (YEDGE[1:] + YEDGE[:-1])/2.
# Bin sizes
XDEL = XEDGE[1]-XEDGE[0]
YDEL = YEDGE[1]-YEDGE[0]

# Pixel area in deg^2; pixel scale is 0.035"/pixel, 
PIXAREA  =((0.035/3600.)*XDEL)*( (0.035/3600.)*YDEL)

# Pre-calculate these instead of doing it in each evaluation of the model...
XX,YY= np.meshgrid(XCENT,YCENT,indexing='ij')
RA,DEC = WCS.wcs_pix2world(XX,YY,1)

# Mask (this is complicated because we need to flatten the arrays)
MASK = IMAGE.contains_points(np.vstack([XX.flatten(),YY.flatten()]).T).T.reshape(XX.shape)
IDX = np.where(MASK==1) # Array indexes

# Initial Eri II model parameters
ERI_NSTAR = 21440  # Richness
ERI_LON = 56.0888  # RA (deg)
ERI_LAT = -43.5304 # Dec (deg)
ERI_EXT = 2.31/60. # Extension (deg)
ERI_ELL = 0.48     # Ellipticity
ERI_PA = 72.6      # Position angle (deg)
ERI_KW = dict(lon=ERI_LON,lat=ERI_LAT,ellipticity=ERI_ELL,
              position_angle=ERI_PA,extension=ERI_EXT)

# Values in pixel coordinates
ERI_X0,ERI_Y0 = WCS.wcs_world2pix(ERI_LON,ERI_LAT,1)
ERI_EXT_PIX = ERI_EXT * 3600 / 0.035 # pixels

def median_interval(data, alpha=0.32):
    """
    Median including Bayesian credible interval.

    Parameters
    ----------
    data  : posterior samples
    alpha : 1 - confidence interval

    Returns
    -------
    [med,[lo, hi]] : median, lower, and upper percentiles
    
    """
    q = [100*alpha/2., 50, 100*(1-alpha/2.)]
    lo,med,hi = np.percentile(data,q)
    return [med,[lo,hi]]

def data(x,y):
    """ Calculate the binned data counts. This only needs to be done
    once (not at each model evaluation), but this seemed easier to
    understand if it paralleled the model counts calculation.

    Parameters
    ----------
    x : the x coordinate of the data
    y : the y coordinate of the data

    Returns
    -------
    data_counts : the data counts in each bin
    """
    data_counts,_,_ = np.histogram2d(x,y,bins=[XEDGE,YEDGE])
    return data_counts

def new_kernel(x,y,lon=ERI_X0,lat=ERI_Y0,ext=ERI_EXT_PIX,ell=ERI_ELL,pa=ERI_PA):
    """ Evaluate the elliptical exponential kernel at coordinates x,y. 
    Normalized to unity over all space...

    Parameters
    ----------
    x: x-coord for evaluating kernel [pix]
    y: y-coord for evaluating kernel [pix]
    lon: x-coord of kernel centroid [pix]
    lat: y-coord of kernel centroid [pix]
    ext: extension [pix]
    ell: ellipticity
    pa:  position angle [deg]

    Returns
    -------
    pdf : probability density (should integrate to unity over all space)
    """

    # Elliptical radius of each x,y coord
    costh = np.cos(np.radians(-pa))
    sinth = np.sin(np.radians(-pa))
    dx = x-lon
    dy = y-lat
    radius = np.sqrt(((dx*costh-dy*sinth)/(1-ell))**2 + (dx*sinth+dy*costh)**2)

    # Exponential radius (re = rh/1.68)
    r_e = ext/1.68 
    #Normalization (integrates to unity over all space?) [stars/pix^2)
    norm = 1./(2*np.pi*r_e**2 * (1-ell) )

    # Exponential PDF
    pdf = norm * np.exp(-radius/r_e)

    return pdf

def model(theta):
    """ Calculate the binned model counts. This extends over the
    entire pixel range, but we will apply the mask later.

    Parameters
    ----------
    theta : the model parameters
    
    Returns
    -------
    model_counts : the model counts in each bin
    """
    richness = theta[0]
    kwargs = dict(lon=theta[1],lat=theta[2])
    # Default values for the other parameters
    kwargs.update(ext=ERI_EXT_PIX,ell=ERI_ELL,pa=ERI_PA)

    # The new kernel in pixel coordinates
    pdf = new_kernel(XX,YY,**kwargs)

    # Calculate the model predicted counts in each pixel
    pixarea = XDEL*YDEL
    model_counts = richness * pdf * pixarea
    return model_counts

def lnlike(theta, x, y):
    """
    Parameters
    ----------
    theta : model parameter array (richness,lon,lat,ext,ell,pa)
    x: x-coordinate of data
    y: y-coordinate of data
    
    Returns
    -------
    lnlike: log-likelihood
    """
    # Calculate the data counts and model predicted counts in each pixel bin
    data_counts = data(x,y)
    model_counts = model(theta)

    # Apply the mask to the data and model. This selects only pixels
    # in the image for calculating the likelihood.
    data_counts_masked = data_counts[IDX]
    model_counts_masked = model_counts[IDX]
    
    # Evaluate Equation C2 from 1912.03302 (ignore k! term)
    lnlike = np.sum(-model_counts_masked + data_counts_masked * np.log(model_counts_masked))
    return lnlike

def lnprior(theta):
    """ The log-prior. Add whatever you want here... 
    
    Parameters
    ----------
    theta : model parameters

    Returns
    -------
    lnprior : log-prior
    """
    rich,x,y = theta[0],theta[1],theta[2]
    if not (10 < rich < 1e5):  return np.inf
    if not (XMIN < x  < XMAX): return np.inf
    if not (YMIN < y  < YMAX): return np.inf
    return 0

def lnprob(theta, x, y):
    """ The log-probability = lnlike + lnprob 

    Parameters
    ----------
    theta : the model parameter vector
    x     : x-coord of the data
    y     : y-coord of the data
    
    Returns
    -------
    lnprob : log-probability
    """
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y)

# Make some diagnostic plots
do_plot=False

# Set to best params from above (same for all sims)
THETA = [ERI_NSTAR,ERI_X0,ERI_Y0]

results = []
hdulist = fits.open('eriII_sims-v3.fits.gz')
for idx,hdu in enumerate(hdulist[1:]):
    print("Reading data and initializing model parameters...")
    cat = hdu.data

    # Remove stars outside the image
    #cat = cat[ (cat['FLAG'] & IMG_BIT) == 0] 

    if do_plot:
        # Example of the masked counts and data
        # (transpose due to difference between imshow and histogram2d...)
        data_counts = data(cat['X'],cat['Y'])
        data_counts_masked = np.copy(data_counts)
        data_counts_masked[np.where(MASK==0)] = np.nan
        plt.imshow(data_counts_masked.T,origin='lower')
        plt.savefig('data_counts_masked.png')
         
        model_counts = model(THETA)
        model_counts_masked = np.copy(model_counts)
        model_counts_masked[np.where(MASK==0)] = np.nan
        plt.imshow(model_counts_masked.T,origin='lower')
        plt.savefig('model_counts_masked.png')
     
    # Initialize and run the mcmc
    print("Running mcmc...")
    ndim, nwalkers = len(THETA), 100
    nthreads,nsamples = 16, 500
    nburn = 100
    pos = [THETA + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
     
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(cat['X'],cat['Y']),
                                    threads=nthreads)
    sampler.run_mcmc(pos,nsamples)
     
    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
     
    rich,[rich_min,rich_max] = median_interval(samples[:,0])
    x,[xmin,xmax] = median_interval(samples[:,1])
    y,[ymin,ymax] = median_interval(samples[:,2])

    # True centroid
    x0,y0 = WCS.wcs_world2pix(hdu.header['LON'],hdu.header['LAT'],1)

    # Results
    res = [x0,y0,x,xmin,xmax,y,ymin,ymax]
    results.append(res)

    if do_plot:
        fig = corner.corner(samples, labels=["rich", "x", "y"])
        fig.savefig("triangle_eri2.png")

results = np.rec.array(results,names=['lon','lat','x','xmin','xmax','y','ymin','ymax'])

filename='results2_b%i_s%i.csv'%(NBINS,nsamples)
print("Writing %s ..."%filename)
rec2csv(results,filename)

