""" Nov 10 wrapper, with r6-r4 flag """

import numpy as np

# from runPhotUtilsDRC_6m4 import runPhotUtils  # run Nov 10
# In the above .py file, there's also a function to create the file structure
# For these files, I have the following directory convention, with the date
# 23Oct. If I decide to re-run things, I will need to change that date.
# Might be good to make a larger intake function in the above so I can
# incorporate that in. Will do that on the next version, as I don't want
# to mess up what I have working so far.
# The directory names (save_dir) modeled after:
#   magCatDir = workDir + '/' + 'drcPhot' + date + '/'
#   catDir = magCatDir + 'catDir_' + targname + '/'
# from getRefDRCfilt import getRefDRCfilt  # run Nov 10
# from drcFiltLinTrans import drcFiltLinTrans  # run Nov 10
# from matchDRCfilt import matchFiltDRC  # run Nov 10
# from makeDRCcmd import makeCMD
# from getJdan import getJdan
# from runPUflcPix import runPhotUtils
# from matchFLCs_first import matchFLCs
# from makeSTDcuts import makeSTDcuts
# from refLinTransFLCfilt import transFLCfilt
# from matchFiltFLC import matchFiltFLC
# from refLinTransFLCdrc import transFLCdrc
# from matchFLCdrc import matchFLCdrc
# from sig2noise import s2n
from plotErrMag31May import plotErr
# from make2CMDinstr import make2CMD
# from outEl import outputEl

# IMPORTANT VARIABLES FOR THIS PROCESS (that can be changed)
work_dir = '../'  # this goes to photRun0520 for me
drcDir = '/Volumes/Spare Data/Hannah_Data/origDRCs/'
jdanDir = '/Volumes/Spare Data/Hannah_Data/'
targFile = work_dir + 'targnamesDirections2.txt'
dateDef = '10Nov'
# drcInfoFile = '/Volumes/Spare Data/Hannah_Data/' + "drcTargInfo_new3.dat"
suffix_ = '_pu.dat'
radius_ = int(4)
matchtol_ = 1  # using for reference star finding
matchtol_f = 2.5  # using for final matching
########

# targname_arr = np.genfromtxt(targFile,dtype='str')
# targname_arr = ['HOROLOGIUM-I']
targname_arr = ['HOROLOGIUM-I','PHOENIX-II','TRIANGULUM-II','RETICULUM-II',
                'SAGITTARIUS-II','HYDRA-II','TUCANA-II']
filt_arr = ['F814W','F606W']

# Run this once, then comment out. May need to cut drcInfoFile
# into smaller pieces to avoid memory overload.
# runPhotUtils(drcInfoFile,radius=radius_,suffix=suffix_,date=dateDef)
# Run the above on 10 Nov @ 11:50 PM

for c1,targname in enumerate(targname_arr):
    # save_dir = work_dir + 'drcPhot' + dateDef + '/' + 'catDir_' \
    #     + targname + '/'
    # save_real = work_dir+ 'drcPhot' + dateDef + '/instCMDs/'

    # getRefDRCfilt(targname,dir=save_dir,matchtol=matchtol_,suffix=suffix_)
    # Run on Nov 10 12:00 AM
    # drcFiltLinTrans(targname,dir=save_dir,suffix=suffix_)
    # linear6d threw an error, but it looks like everything was completed.
    # There's now an extra drcFiltTrans.dat in directories from a past run.
    # Run on Nov 10 12:00 AM
    # matchFiltDRC(targname,dir=save_dir,matchtol=matchtol_f,suffix=suffix_)
    # Run on Nov 10 12:00 AM
    # makeCMD(targname,dir=save_dir)

    # for c2, filt in enumerate(filt_arr):
        # jdan,drc_file = getJdan(targname,filt,dir=jdanDir,drc=True)
        # runPhotUtils(targname,filt,radius=radius_,jdanUse=jdan,
        #              drcFile=drcDir + drc_file,saveDir=save_dir,flcDir='../',
        #              suffix='_pu.dat')  # running on Nov 12 @ 12:20 AM
        # print('Matching FLCs for {0} {1}...'.format(targname,filt))
        # matchFLCs(targname,filt,jdan,dir=save_dir,matchtol=matchtol_f)
        # print('Making STD cuts for {0} {1}...'.format(targname,filt))
        # makeSTDcuts(filt,dir=save_dir)

    # transFLCfilt(targname,dir=save_dir,matchtol=matchtol_)
    # print('Matching FLCs between filters for {0}...'.format(targname))
    # matchFiltFLC(targname,dir=save_dir,matchtol=matchtol_)
    # flc_file = targname + '_allMatchedFLC.dat'
    # drc_file = targname + '_matchedDRCfilt.dat'
    # print('Transforming FLC to DRC for {0}...'.format(targname))
    # transFLCdrc(targname,flcFile=flc_file,drcFile=drc_file,dir=save_dir,
    #             matchtol=matchtol_)
    # flc_file = 'flcDRCtrans_' + targname + '.dat'
    # print('Matching FLCs to DRC for {0}...'.format(targname))
    # matchFLCdrc(targname,flcFile=flc_file,drcFile=drc_file,dir=save_dir,
    #             matchtol=matchtol_)
    # s2n(targname,dir=save_dir)
    # plotErr(targname,dir=save_dir)
    plotErr(targname)
    # make2CMD(targname,dir=save_dir,saveDir=save_real)
    # outputEl(targname,dir=save_dir)
