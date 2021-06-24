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
match_files_dir = '/Volumes/Spare Data/Hannah_Data/drc_flc_match/'
flc_dir = '/Volumes/Spare Data/Hannah_Data/hor1dir6pix/'

# assign column numbers to the columns for the first file
RA_f606w,DEC_f606w,flux_f606w,c_star_f606w,xr_1_f606w,yr_1_f606w,xt_1_f606w,yt_1_f606w,mag1_f606w,xr_2_f606w,yr_2_f606w,xt_2_f606w,yt_2_f606w,mag2_f606w,xr_3_f606w,yr_3_f606w,xt_3_f606w,yt_3_f606w,mag3_f606w,xr_4_f606w,yr_4_f606w,xt_4_f606w,yt_4_f606w,mag4_f606w,mean_f606w,stdev_f606w,pos_std_f606w,cut_flag_f606w,cut_idx_f606w,num_abv_std_f606w,RA_f814w,DEC_f814w,flux_f814w,c_star_f814w,xr_1_f814w,yr_1_f814w,xt_1_f814w,yt_1_f814w,mag1_f814w,xr_2_f814w,yr_2_f814w,xt_2_f814w,yt_2_f814w,mag2_f814w,xr_3_f814w,yr_3_f814w,xt_3_f814w,yt_3_f814w,mag3_f814w,xr_4_f814w,yr_4_f814w,xt_4_f814w,yt_4_f814w,mag4_f814w,mean_f814w,stdev_f814w,pos_std_f814w,cut_flag_f814w,cut_idx_f814w,num_abv_std_f814w = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59


# assign column numbers for the second file (in this case the master file)
RA_v,DEC_v,x_v,y_v,fAper_v,fErr_v,magAper_v,magErr_v,magRaw_v,magRed_v,magAbs_v,elong_v,ellip_v,class_Star_v,RA_i,DEC_i,x_i,y_i,fAper_i,fErr_i,magAper_i,magErr_i,magRaw_i,magRed_i,magAbs_i,elong_i,ellip_i,class_Star_i,corrF_errV,corrF_errI,corrM_errV,corrM_errI = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31

# List of sources already matched to the master frame with x/y in the original frame
match = np.loadtxt(match_files_dir+"flc_upper_hor1_1e4.dat", comments="#")

# List of matched sources in the master frame
master = np.loadtxt(match_files_dir+"drc_upper_hor1_1e4.dat", comments="#")

# List of all sources in original frame for the transformation to be applied to
all = np.loadtxt(flc_dir+"HORI_pix_2212_d6_dist.dat", comments="#")

# Generate an equal length array containing only 1 for the weights
weights = np.zeros((len(match)))
weights.fill(1.0)

# Call the linear transformation code with the following order:
# x-pos original, y-pos original, x-pos master, y-pos master,
# x-pos weights, y-pos weights, x-pos all sources, y-pos all sources

# new_match, new_all = test_linear(match[:,xt_1_f606w], match[:,yt_1_f606w], master[:,x_v], master[:,y_v], weights, weights, all[:,xt_1_f606w], all[:,yt_1_f606w])

new_match, new_all = test_linear(match[:,xt_1_f814w], match[:,yt_1_f814w], master[:,x_v], master[:,y_v], weights, weights, all[:,xt_1_f814w], all[:,yt_1_f814w])

# Save out the two new transformed x/y positions
np.savetxt(match_files_dir+"upper_hor1_814_xy.dat", new_all, fmt="%1.6f")
