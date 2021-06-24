""" Hopefully a cleaner, more straightforward wrapper """

import numpy as np

# from runPhotUtilsDRC import runPhotUtils  # run Oct 26
# In the above .py file, there's also a function to create the file structure
# For these files, I have the following directory convention, with the date
# 23Oct. If I decide to re-run things, I will need to change that date.
# Might be good to make a larger intake function in the above so I can
# incorporate that in. Will do that on the next version, as I don't want
# to mess up what I have working so far.
# The directory names (save_dir) modeled after:
#   magCatDir = workDir + '/' + 'drcPhot' + date + '/'
#   catDir = magCatDir + 'catDir_' + targname + '/'
# from getRefDRCfilt import getRefDRCfilt  # run Oct 26
# from drcFiltLinTrans import drcFiltLinTrans  # run Oct 26
# from matchDRCfilt import matchFiltDRC  # run Oct 27
from getJdan import getJdan


# IMPORTANT VARIABLES FOR THIS PROCESS (that can be changed)
work_dir = '../'  # this goes to photRun0520 for me
drcDir = '/Volumes/Spare Data/Hannah_Data/origDRCs/'
jdanDir = '/Volumes/Spare Data/Hannah_Data/'
targFile = work_dir + 'targnamesDirections2.txt'
dateDef = '23Oct'
drcInfoFile = '/Volumes/Spare Data/Hannah_Data/' + "drcTargInfo_new3.dat"
suffix_ = '_pu.dat'
radius_ = int(4)
matchtol_ = 1  # using for reference star finding
matchtol_f = 2.5  # using for final matching
########

targname_arr = np.genfromtxt(targFile,dtype='str')
filt_arr = ['F814W','F606W']

# Run this once, then comment out. May need to cut drcInfoFile
# into smaller pieces to avoid memory overload.
# runPhotUtils(drcInfoFile,radius=radius_,suffix=suffix_)
# Run the above on October 22

for c1,targname in enumerate(targname_arr):
    save_dir = work_dir + '/' + 'drcPhot' + dateDef + '/' + 'catDir_' \
        + targname + '/'

    # getRefDRCfilt(targname,dir=save_dir,matchtol=matchtol_,suffix=suffix_)
    # Run on Oct 26 @ 2 PM
    # drcFiltLinTrans(targname,dir=save_dir,suffix=suffix_)
    # linear6d threw an error, but it looks like everything was completed.
    # There's now an extra drcFiltTrans.dat in directories from a past run.
    # Run on Oct 26 @ 8 PM
    # matchFiltDRC(targname,dir=save_dir,matchtol=matchtol_f,suffix=suffix_)
    # Run on Oct 27 @ 12 AM

    for c2, filt in enumerate(filt_arr):
        jdan = getJdan(targname,filt,dir=jdanDir)

#
#
#
