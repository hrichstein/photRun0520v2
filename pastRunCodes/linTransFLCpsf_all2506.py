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

c_star_f606w, mag1_f606w, mag2_f606w, mag3_f606w, mag4_f606w, xr1_f606w, yr1_f606w, xt1_f606w, yt1_f606w = 0, 1, 2, 3, 4, 5, 6, 7, 8


RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

# F606W
# flcCat = np.genfromtxt(dir+'matched_w_MagsPos1106r2.dat')
# drcCat = np.genfromtxt(dir +'drc_useful_f606w.dat')
#
# flcArr = np.array([flcCat[484],flcCat[184],flcCat[3],flcCat[2069],flcCat[150],flcCat[1097],flcCat[804],flcCat[789],flcCat[63]])
#
# drcArr = np.array([drcCat[338],drcCat[115],drcCat[0],drcCat[1282],drcCat[47],drcCat[834],drcCat[581],drcCat[464],drcCat[50]])

flcCat = np.genfromtxt(dir+'matchedFLCdrc2506_comb.dat')
psfCat = np.genfromtxt(dir +'matchedPSFaper1706_tc.dat')

xFLC = np.array([1942.41,2733.6358,3210.4213,275.484,1854.0268,3119.8681,3321.6462,1012.2739,2960.6373,1586.7523])

yFLC = np.array([336.0805,2591.7028,721.6712,2016.7489,2635.92,3469.8079,1810.6803,290.3904,3657.488,1715.3336])

xPSF = np.array([2434.9197,4718.1953,2887.3555,4013.6179,4714.2852,5610.8672,3974.3972,2339.0757,5788.7881,3783.894])

yPSF = np.array([5232.4946,4569.6943,3994.6946,6979.1934,5445.4702,4233.6597,3943.2954,6153.9668,4402.1094,5660.8569])


match_arr = np.zeros((len(xFLC),2))

match_arr[:,0] = xFLC
match_arr[:,1] = yFLC

master_arr = np.zeros((len(xPSF),2))

master_arr[:,0] = xPSF
master_arr[:,1] = yPSF


def linTrans(targname,filt,dir='catRawMags1305/catDir/'):
# matched_w_MagsPos2705.dat
    # match = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    # master = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    match = match_arr
    master= master_arr

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    all = flcCat

    # outname = 'flcDRCmatch2506_f606w.dat'
    outname = 'flcPSFmatch2506_all.dat'

    new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xt1_f606w],all[:,yt1_f606w])

    np.savetxt(dir+outname, new_all, fmt="%1.6f")


    return None

def addTranscols(targname,filt,dir='catRawMags1305/catDir/'):

    # cat = np.genfromtxt(dir+'matched_w_MagsPos1106r2.dat')
    # transCat = np.genfromtxt(dir + 'flcDRCmatch2506_f606w.dat')

    cat = np.genfromtxt(dir+'matchedFLCdrc2506_comb.dat')
    transCat = np.genfromtxt(dir + 'flcPSFmatch2506_all.dat')


    newCol = np.zeros((len(cat),2))

    newCol[:,0] = transCat[:,0]
    newCol[:,1] = transCat[:,1]

    cat = np.hstack((cat, newCol))

    header = 'c_star_f606w mag1_f606w mag2_f606w mag3_f606w mag4_f606w xr1_f606w yr1_f606w xt1_f606w yt1_f606w xDRC_trans_f606w yDRC_trans_f606w xDRC_mat_f606w yDRC_mat_f606w magDRC_f606w id_cat_f606w mean_f606w stdev_f606w cut_flag_f606w idx_cut_f606w num_abv_std_f606w magZPT_f606w magZPTerr_f606w c_star_f814w mag1_f814w mag2_f814w mag3_f814w mag4_f814w xr1_f814w yr1_f814w xt1_f814w yt1_f814w xDRC_trans_f814w yDRC_trans_f814w xDRC_mat_f814w yDRC_mat_f814w magDRC_f814w id_cat_f814w mean_f814w stdev_f814w cut_flag_f814w idx_cut_f814w num_abv_std_f814w magZPT_f814w magZPTerr_f814w xPSF yPSF'

    # np.savetxt(dir+'flcDRCpos2506_f606w.dat', cat, header=header)

    np.savetxt(dir+'flcPSFpos2506_all.dat', cat, header=header)
    return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'

linTrans(targname,filt,dir='catRawMags1305/catDir/')
addTranscols(targname,filt,dir='catRawMags1305/catDir/')


# Iter      18    CHI-SQUARE =  3.710903345  DOF =  12
#    P0 = -21.88823476
#    P1 = 0.9935535276
#    P2 = 0.03846484148
#    P3 = 156.8124343
#    P4 = -0.0388812208
#    P5 = 0.9934957381
