'''

hst_fun.py

This is a local library file containing functions pertaining to
manipulating HST data.

'''

import numpy as np
from astropy.io import fits

###########################################################################
### This is the function to perform the distortion correction, adapted
### from Tony Sohn's IDL implementation of Jay Anderson's Fortran
### distortion correction. It takes three input arguments: "filebase"
### which is the location of the distortion coeficients (stopping short
### of the gcx.fits/gcy.fits; "xr" which is a numpy array (of shape (n,))
### with the x-coordinates of the sources; "yr" which is the same as "xr"
### but for the y-coordinates.
### This function also assumes that the image being corrected
### is a 4096 x 4096 image.

def wfc3uvis(filebase, xr, yr):

### Create two new arrays of integers to interact with the solution files
  ix = np.zeros(len(xr))
  iy = np.zeros(len(yr))

  for i in range(len(xr)):
    ix[i] = int(xr[i])
    iy[i] = int(yr[i])


### Load the distortion solution coefficients
  gcx = fits.open(filebase+"_gcx.fits")
  gcxim = gcx[0].data
  gcx.close()

  gcy = fits.open(filebase+"_gcy.fits")
  gcyim = gcy[0].data
  gcy.close()


### Check for oddities in the associated source indicies

  for i in range(len(ix)):
    if (ix[i] < 1):
      ix[i] = 1
    elif (ix[i] > 4095):
      ix[i] = 4095

  for i in range(len(iy)):
    if (iy[i] < 1):
      iy[i] = 1
    elif (iy[i] > 4095):
      iy[i] = 4095


### Create differentials for the correction
  fx = xr - ix
  fy = yr - iy

### Create array to store the corrected values
  cor = np.zeros((len(xr),2))

### Perform the distortion correction

  for i in range(len(xr)):
    cor[i][0] = (1-fx[i])*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcxim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcxim[int(iy[i])][int(ix[i])]

    cor[i][1] = (1-fx[i])*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcyim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcyim[int(iy[i])][int(ix[i])]

###Return the new arrays
  return cor

###########################################################################
def acsDistortion(filebase, xr, yr):

### Create two new arrays of integers to interact with the solution files
  ix = np.zeros(len(xr))
  iy = np.zeros(len(yr))

  for i in range(len(xr)):
    ix[i] = int(xr[i])
    iy[i] = int(yr[i])


### Load the distortion solution coefficients
  gcx = fits.open(filebase+"_gcx_SM4.fits")
  gcxim = gcx[0].data
  gcx.close()

  gcy = fits.open(filebase+"_gcy_SM4.fits")
  gcyim = gcy[0].data
  gcy.close()


### Check for oddities in the associated source indicies

  for i in range(len(ix)):
    if (ix[i] < 1):
      ix[i] = 1
    elif (ix[i] > 4095):
      ix[i] = 4095

  for i in range(len(iy)):
    if (iy[i] < 1):
      iy[i] = 1
    elif (iy[i] > 4095):
      iy[i] = 4095


### Create differentials for the correction
  fx = xr - ix
  fy = yr - iy

### Create array to store the corrected values
  cor = np.zeros((len(xr),2))

### Perform the distortion correction

  for i in range(len(xr)):
    cor[i][0] = (1-fx[i])*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcxim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcxim[int(iy[i])][int(ix[i])]

    cor[i][1] = (1-fx[i])*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcyim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcyim[int(iy[i])][int(ix[i])]

###Return the new arrays
  return cor



def get_gcd(ra1, dec1, ra2, dec2):
    """
    (Private internal function)
    Returns great circle distance.  Inputs in degrees.

    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.
    """
    from numpy import radians, degrees, sin, cos, arctan2, hypot

    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"
    lambs = radians(ra1)
    phis = radians(dec1)
    lambf = radians(ra2)
    phif = radians(dec2)

    dlamb = lambf - lambs

    numera = cos(phif) * sin(dlamb)
    numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
    numer = hypot(numera, numerb)
    denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)


    return degrees(arctan2(numer, denom))





    #End
