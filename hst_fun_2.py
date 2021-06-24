'''

hst_fun_2.py

This is a local library file containing functions pertaining to
manipulating HST data.

'''

import numpy as np
from astropy.io import fits

###########################################################################
### This is the function to perform the distortion correction for HST
### magnitudes, adapted from Tony Sohn's IDL implementation of
### Jay Anderson's Fortran distortion correction for the spatial locations
### of stars. It takes three input arguments: "filebase"
### which is the location of the distortion coeficients (stopping short
### of the gcx.fits/gcy.fits; "xr" which is a numpy array (of shape (n,))
### with the x-coordinates of the sources; "yr" which is the same as "xr"
### but for the y-coordinates.
### This function also assumes that the image being corrected
### is a 4096 x 4096 image.
###
### It returns a list of corrections to be applied to the
### instrumental magnitudes based on the location of the sources.

def hst_gc_z(filebase, xr, yr):

### Create two new arrays of integers to interact with the solution files
  ix = np.zeros(len(xr))
  iy = np.zeros(len(yr))

  for i in range(len(xr)):
    ix[i] = int(xr[i])
    iy[i] = int(yr[i])


### Load the distortion solution coefficients
  gcz = fits.open(filebase+"_gcz_SM4.fits")
  gczim = gcz[0].data
  gcz.close()

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
  cor = np.zeros((len(xr),1))

  iy = np.asarray(iy,dtype=int)
  ix = np.asarray(ix,dtype=int)

  fx = np.asarray(fx,dtype=int)
  fy = np.asarray(fy,dtype=int)
### Perform the distortion correction

  for i in range(len(xr)):
    cor[i][0] = (1-fx[i])*(1-fy[i])*gczim[int(iy[i])-1][int(ix[i])-1]+ \
       (1-fx[i])*( fy[i] )*gczim[int(iy[i])][int(ix[i])-1] + \
       ( fx[i] )*(1-fy[i])*gczim[int(iy[i])-1][int(ix[i])] + \
       ( fx[i] )*( fy[i] )*gczim[int(iy[i])][int(ix[i])]

    # print(cor)
  # for i in range(len(xr)):
  #   cor[i][0] = (1-fx[i])*(1-fy[i])*gczim[(iy[i])-1][(ix[i])-1] + \
  #          (1-fx[i])*( fy[i] )*gczim[(iy[i])][(ix[i])-1] + \
  #          ( fx[i] )*(1-fy[i])*gczim[(iy[i])-1][(ix[i])] + \
  #          ( fx[i] )*( fy[i] )*gczim[(iy[i])][(ix[i])]
  #
  # for i in range(4):
  #
  #     print(int(iy[i]))
  #     print(int(ix[i]))

###Return the new arrays
  return cor
  # return None
###########################################################################
