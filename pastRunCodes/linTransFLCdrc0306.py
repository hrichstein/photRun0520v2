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

# xt, yt = 11,12 # column numbers in dc file

RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

# Original FLC X,Y
x_arr = np.array([1905.8447, 2733.5304, 1735.3688, 1586.7118, 1870.5514, 1578.4376,1843.0081, 1942.4414])
y_arr = np.array([3323.9234, 2591.8753, 1915.1304,  1715.1404, 1625.8418,  736.2737, 625.12, 336.7415])

# Rotated FLC X,T
# x_arr = np.array([3489.19255567, 2838.66038074, 2070.59434489, 1857.44664442,1795.37820268,  882.18082431,  796.53370475,  518.8455224 ])
# y_arr= np.array([-1583.11278784, -2476.29022154, -1546.57012504,
#  -1417.4811784,  -1708.4909154,  -1501.77268973, -1775.66548884,
#  -1901.9129197])

# Rotated FLC X,T with 1899, 7065 offset added; resulting from first run matching Rotated FLC X,T with 1899, 7065 offset added to PSF
# x_arr = np.array([5484.295887,   4802.40091105, 4076.11897452, 3869.44630487, 3796.21169842,2896.99915618, 2801.01482972, 2520.06805614])
# y_arr = np.array([5543.16195301, 4681.00524382, 5635.32278369, 5772.019284,   5485.21671074,5726.62768698, 5457.76641413, 5343.21476812])

# Rotated FLC X,T, no offset applied; resulting from first run matching Rotated FLC X,T with 1899, 7065 offset added to PSF
# x_arr = np.array([3585.295887,   2903.40091105, 2177.11897452, 1970.44630487, 1897.21169842,997.99915618,  902.01482972,  621.06805614])
# y_arr = np.array([-1521.83804699, -2383.99475618, -1429.67721631, -1292.980716,
#  -1579.78328926, -1338.37231302, -1607.23358587, -1721.78523188])

match_arr = np.zeros((len(x_arr),2))

match_arr[:,0] = x_arr
match_arr[:,1] = y_arr


# PSF positions
x_mas = np.array([5400.123, 4718.1953, 3991.6089, 3783.894, 3711.5825, 2812.5161,2716.2554, 2434.9197])

y_mas = np.array([5431.417,  4569.6943, 5523.9956, 5660.8569, 5374.0845, 5615.8149,5347.2842, 5232.4946])

master_arr = np.zeros((len(x_mas),2))

master_arr[:,0] = x_mas
master_arr[:,1] = y_mas

def linTrans(targname,filt,dir='catRawMags1305/catDir/'):
# matched_w_MagsPos2705.dat
    # match = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    # master = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    match = match_arr
    master= master_arr

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    jdanUse = getJdan(targname,filt)
    #
    # for ff in range(len(jdanUse)):
    #
    #     if ff==0:
    #         xr_x = xt1
    #         yr_x = yt1
    #
    #     elif ff==1:
    #         xr_x = xt2
    #         yr_x = yt2
    #
    #     elif ff==2:
    #         xr_x = xt3
    #         yr_x = yt3
    #
    #     else:
    #         xr_x = xt4
    #         yr_x = yt4

        # all = np.loadtxt(dir+"master_ids_"+targname+"_"+filt+"_"+str(iter-1)+".dat")
    all = np.loadtxt(dir+"matched_w_MagsPos2705r2.dat")
    # all = np.loadtxt(dir+'flcRo2PSF_round2.dat')
    # all = np.loadtxt(dir+'jdan21l8q_HOROLOGIUM-I_F814W_3105DRC2PSFall.dat')

    # all1000_idx = np.argsort(all[:,7])
    # #
    # all = all[all1000_idx]

    outname = "flc2PSF_round1.dat"
    out2 = "flc2PSF_round1_newMatch.dat"

    # new_match, new_all = test_linear(match[:,0], match[:,1], master[:,0], master[:,1], weights, weights, all[:,xt], all[:,yt])

    # new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xt1],all[:,yt1])

    new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xt1],all[:,yt1])

    np.savetxt(dir+outname, new_all, fmt="%1.6f")
    np.savetxt(dir+out2, new_match, fmt="%1.6f")


    return None

# def addTranscols(targname,filt,dir='catRawMags1305/catDir/'):
#
#     jdanUse = getJdan(targname,filt)
#
#     for ff in range(len(jdanUse)):
#
#         cat = np.loadtxt(dir+jdanUse[ff]+"_"+targname+"_"+filt+"_dc.dat")
#
#         transCat = np.loadtxt(dir+jdanUse[ff]+"_"+targname+"_"+filt+"_t_2705r3.dat")
#
#         newCol = np.zeros((len(cat),2))
#
#         newCol[:,0] = transCat[:,0]
#         newCol[:,1] = transCat[:,1]
#
#         cat = np.hstack((cat, newCol))
#
#         header = 'flags RA DEC xr yr flux c_star magr id xc yc xt yt'
#
#         np.savetxt(dir+jdanUse[ff] + "_"+targname+"_"+filt+"_at_2705rdrc.dat", cat, fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.6f %1.6f", header=header)
#
#
#     return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'

linTrans(targname,filt,dir='catRawMags1305/catDir/')
# addTranscols(targname,filt,dir='catRawMags1305/catDir/')
