import numpy as np
from linear6d import *
from getJdan import *

dir = 'catRawMags1305/catDir/'
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

# xt,yt in this case are just the offset corrected, distortion corrected positions

xt, yt = 11,12 # column numbers in dc file

RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

def linTrans(targname,filt,dir='catRawMags1305/catDir/'):
# matched_w_MagsPos2705.dat
    match = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    master = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):

        if ff==0:
            xr_x = xt1
            yr_x = yt1

        elif ff==1:
            xr_x = xt2
            yr_x = yt2

        elif ff==2:
            xr_x = xt3
            yr_x = yt3

        else:
            xr_x = xt4
            yr_x = yt4

        # all = np.loadtxt(dir+"master_ids_"+targname+"_"+filt+"_"+str(iter-1)+".dat")
        all = np.loadtxt(dir+jdanUse[ff]+"_"+targname+"_"+filt+"_at_2705r2.dat")

        outname = jdanUse[ff]+"_"+targname+"_"+filt+"_t_2705r3.dat"

        new_match, new_all = test_linear(match[:,xr_x], match[:,yr_x], master[:,xt1], master[:,yt1], weights, weights, all[:,xt], all[:,yt])

        np.savetxt(dir+outname, new_all, fmt="%1.6f")


    return None

def addTranscols(targname,filt,dir='catRawMags1305/catDir/'):

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):

        cat = np.loadtxt(dir+jdanUse[ff]+"_"+targname+"_"+filt+"_dc.dat")

        transCat = np.loadtxt(dir+jdanUse[ff]+"_"+targname+"_"+filt+"_t_2705r3.dat")

        newCol = np.zeros((len(cat),2))

        newCol[:,0] = transCat[:,0]
        newCol[:,1] = transCat[:,1]

        cat = np.hstack((cat, newCol))

        header = 'flags RA DEC xr yr flux c_star magr id xc yc xt yt'

        np.savetxt(dir+jdanUse[ff] + "_"+targname+"_"+filt+"_at_2705r3.dat", cat, fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.6f %1.6f", header=header)


    return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'

linTrans(targname,filt,dir='catRawMags1305/catDir/')
addTranscols(targname,filt,dir='catRawMags1305/catDir/')
