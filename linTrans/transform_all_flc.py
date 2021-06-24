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

'''

import numpy as np
from linear6d import *

####

# Main directory for the project I'm working on
dir = "/Volumes/Spare Data/Hannah_Data/Linear_transform_code/"
main_dir = '/Volumes/Spare Data/Hannah_Data/hor1dir2804/'

def trans_all_flc(filter):
# List of sources already matched to the master frame with x/y in the original frame

    match = np.genfromtxt(main_dir+'flcs_w_drcMAGS_'+filter+'.dat',names=True)

# List of matched sources in the master frame
    master = np.genfromtxt(main_dir+'flcs_w_drcMAGS_'+filter+'.dat',names=True)

# List of all sources in original frame for the transformation to be applied to
    all = np.genfromtxt(main_dir+"HORI_comb2804.dat", names=True)

# Generate an equal length array containing only 1 for the weights
    weights = np.zeros((len(match)))
    weights.fill(1.0)

    # Call the linear transformation code with the following order:
    # x-pos original, y-pos original, x-pos master, y-pos master,
    # x-pos weights, y-pos weights, x-pos all sources, y-pos all sources

    new_match, new_all = test_linear(match["xt_1_"+filter], match["yt_1_"+filter], master["xdrc_up_"+filter], master["ydrc_up_"+filter], weights, weights, all["xt_1_"+filter], all["yt_1_"+filter])

    # Save out the two new transformed x/y positions
    np.savetxt(main_dir+"lit_all_flc_up_"+filter+"_xy.dat", new_all, fmt="%1.6f")

    return None

trans_all_flc('f814w')
