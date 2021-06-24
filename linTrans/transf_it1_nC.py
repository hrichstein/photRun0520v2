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

import numpy as np
from linear6d import *

workdir = "/Volumes/Spare Data/Hannah_Data/"

def getJdan(targname,filter):

    jdanList_file = np.loadtxt(workdir+targname+"_flcs.txt",dtype="<U50")

    jdan = jdanList_file[:,0]
    jdanFilter = jdanList_file[:,1]

    jdanUse = []
    for jj in range(len(jdan)):
        if jdanFilter[jj]==filter:
            jdanUse.append(jdan[jj])

    jdanUse = np.array(jdanUse)

    return jdanUse

####
# Main directory for the project I'm working on
dir = "/Volumes/Spare Data/Hannah_Data/Linear_transform_code/"
# match_files_dir = '/Volumes/Spare Data/Hannah_Data/drc_flc_match/'
# flc_dir = '/Volumes/Spare Data/Hannah_Data/hor1dir6pix/'

match_files_dir = '/Volumes/Spare Data/Hannah_Data/hor1dir3103/'
flc_dir = '/Volumes/Spare Data/Hannah_Data/hor1flcs/catRawMags1703/catDir/'

# assign column numbers to the columns for the file with sources to be transformed
xr, yr, flux, c_star, magr, id, xc, yc = 0, 1, 2, 3, 4, 5, 6, 7


# assign column numbers for the second file (in this case the master file)

RA, DEC, flux, c_star, xr_1, yr_1, xc_1, yc_1, xt_1, yt_1, mag1, xr_2, yr_2, xc_2, yc_2, xt_2, yt_2, mag2, xr_3, yr_3, xc_3, yc_3, xt_3, yt_3, mag3, xr_4, yr_4, xc_4, yc_4, xt_4, yt_4, mag4, mean, stdev, pos_std, cut_flag, cut_idx, num_abv_std = 0, 1, 2, 3, 4, 5, 6, 7, 8, \
    9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37

# xr, yr, flux, c_star, magr, id, xc, yc = 0, 1, 2, 3, 4, 5, 6, 7
# List of sources already matched to the master frame with x/y in the original frame

def trans_flc2flc(targname,filter):

    match = np.loadtxt(match_files_dir+"HOR-I_refStars_"+filter+"_nC.dat", comments="#")

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    master = np.loadtxt(match_files_dir+"HOR-I_refStars_"+filter+"_nC.dat", comments="#")

    jdanUse = getJdan(targname,filter)

    for ff in range(len(jdanUse)):

        if ff==0:
            xr_x = xc_1
            yr_x = yc_1

        elif ff==1:
            xr_x = xc_2
            yr_x = yc_2

        elif ff==2:
            xr_x = xc_3
            yr_x = yc_3

        else:
            xr_x = xc_4
            yr_x = yc_4

        all = np.loadtxt(flc_dir+jdanUse[ff]+"_"+targname+"_"+filter+"_dc.dat", comments="#")

        outname = jdanUse[ff]+"_"+targname+"_"+filter+"_tnC.dat"

        new_match, new_all = test_linear(match[:,xr_x], match[:,yr_x], master[:,xc_1], master[:,yc_1], weights, weights, all[:,xc], all[:,yc])

        # print('working')
        np.savetxt(match_files_dir+outname, new_all, fmt="%1.6f")

    return None
# The xrs across the board
# List of all sources in original frame for the transformation to be applied to


# Generate an equal length array containing only 1 for the weights
trans_flc2flc('HOROLOGIUM-I','F606W')
trans_flc2flc('HOROLOGIUM-I','F814W')

# Call the linear transformation code with the following order:
# x-pos original, y-pos original, x-pos master, y-pos master,
# x-pos weights, y-pos weights, x-pos all sources, y-pos all sources

# new_match, new_all = test_linear(match[:,xt_1_f606w], match[:,yt_1_f606w], master[:,x_v], master[:,y_v], weights, weights, all[:,xt_1_f606w], all[:,yt_1_f606w])
#


# Save out the two new transformed x/y positions
#
