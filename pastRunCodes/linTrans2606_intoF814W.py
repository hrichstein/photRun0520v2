import numpy as np
from linear6d import *

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

dir = 'sDRC_2606/'

arr = np.genfromtxt(dir+'matchedDRCref2606.dat')

f606w_cat = np.genfromtxt(dir+'cat_HOR-I_F606W_cStarCut.dat')

flags_f606w, RA_f606w, DEC_f606w, xr_f606w, yr_f606w, flux_f606w, c_star_f606w, magr_f606w, id_f606w, flags_f814w, RA_f814w, DEC_f814w, xr_f814w, yr_f814w, flux_f814w, c_star_f814w, magr_f814w, id_f814w = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17

flags, RA, DEC, xr, yr, flux, c_star, magr, id = 0, 1, 2, 3, 4, 5, 6, 7, 8

match_arr = np.zeros((len(arr),2))

match_arr[:,0] = arr[:,xr_f606w]
match_arr[:,1] = arr[:,yr_f606w]

master_arr = np.zeros((len(arr),2))

master_arr[:,0] = arr[:,xr_f814w]
master_arr[:,1] = arr[:,yr_f814w]


def linTrans(targname,filt):
# matched_w_MagsPos2705.dat
    # match = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    # master = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    match = match_arr
    master= master_arr

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    all = f606w_cat

    # outname = 'flcDRCmatch2506_f606w.dat'
    outname = 'drcDRCmatch2606_f814w.dat'

    new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xr],all[:,yr])

    np.savetxt(dir+outname, new_all, fmt="%1.6f")
    # np.savetxt(dir+'flcDRC0406newMatch.dat',new_match,fmt='%1.6f')


    return None

def addTranscols(targname,filt):

    # cat = np.genfromtxt(dir+'matched_w_MagsPos1106r2.dat')
    # transCat = np.genfromtxt(dir + 'flcDRCmatch2506_f606w.dat')

    cat = np.genfromtxt(dir+'cat_HOR-I_F606W_cStarCut.dat')
    transCat = np.genfromtxt(dir + 'drcDRCmatch2606_f814w.dat')

    newCol = np.zeros((len(cat),2))

    newCol[:,0] = transCat[:,0]
    newCol[:,1] = transCat[:,1]

    cat = np.hstack((cat, newCol))

    header = 'flags_f606w RA_f606w DEC_f606w xr_f606w yr_f606w flux_f606w c_star_f606w magr_f606w id_f606w xr_f814w_trans yr_f814w_trans'

    # np.savetxt(dir+'flcDRCpos2506_f606w.dat', cat, header=header)

    np.savetxt(dir+'drcDRCpos2606_intoF814W.dat', cat, header=header)
    return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'

linTrans(targname,filt)
addTranscols(targname,filt)


# Iter      16    CHI-SQUARE =  0.01621569722  DOF =  20
#    P0 = -0.06300029445
#    P1 = 1.000015683
#    P2 = 1.669844601e-05
#    P3 = -0.01987557594
#    P4 = 4.015278338e-06
#    P5 = 1.000005038
