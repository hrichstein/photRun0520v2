import numpy as np
from linear6d import *
from getJdan import *

'''

transform.py

This code serves as the main interface for using the mpfit Python code (mpfit.py) to
determine the coefficients for a 6-parameter linear transformation (linear6d.py).
For the input, it takes 8 different 1D arrays:
- x and y positions of a list of sources to be transformed into another frame that
have already been matched to a corresponding list of sources in the other frame.

- x and y positions of the corresponding list in the intended frame

- weights for the sources (this is currently set to 1.0 for every source
to give equal weight

- x and y positions of all sources that are to be transformed into another frame. Once
the 6-param coefficients have been determined using the above matched sources, the
transformation is applied to this full list of sources and is saved out.


The first set of coordinates is the reference sources in your particular secondary dither. The second set of coordinates is the reference sources in the primary dither. And the third set is all sources in your secondary dither.
'''

# flags, RA, DEC, xr, yr, flux, c_star, magr, id = = 0, 1, 2, 3, 4, 5, 6, 7, 8

xc, yc = 9, 10 # column numbers in dc file

wcsRA, wcsDEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xc1, yc1, xt1, yt1, xr2, yr2, xc2, yc2, xt2, yt2, xr3, yr3, xc3, yc3, xt3, yt3, xr4, yr4, xc4, yc4, xt4, yt4, mean, stdev, pix_std, coord_std, cut_flag, idx_cut, num_abv_std  = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47


def linTrans(targname,filt,iter,workDir='./'):

    match = np.loadtxt(workDir+targname+"_refStars_"+filt+'_'+str(iter-1)+'.dat')
    master = np.loadtxt(workDir+targname+"_refStars_"+filt+'_'+str(iter-1)+'.dat')

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):

        if ff==0:
            xr_x = xc1
            yr_x = yc1

        elif ff==1:
            xr_x = xc2
            yr_x = yc2

        elif ff==2:
            xr_x = xc3
            yr_x = yc3

        else:
            xr_x = xc4
            yr_x = yc4

        # all = np.loadtxt(workDir+"master_ids_"+targname+"_"+filt+"_"+str(iter-1)+".dat")
        all = np.loadtxt(workDir+jdanUse[ff]+"_"+targname+"_"+filt+"_dc.dat")

        outname = jdanUse[ff]+"_"+targname+"_"+filt+"_t_"+str(iter)+".dat".format(iter)

        new_match, new_all = test_linear(match[:,xr_x], match[:,yr_x], master[:,xc1], master[:,yc1], weights, weights, all[:,xc], all[:,yc])

        np.savetxt(workDir+outname, new_all, fmt="%1.6f")


    return None
