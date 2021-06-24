import numpy as np
import os
from astropy.io import fits
import time

from getJdan import getJdan
from exportDRCuseful0707 import *
from f2mag0707 import *
# from hst_func import *
# from linTrans import *
# from f2mag0707 import f2mag_dirs
from initialCorrMatch0707 import *
from linTrans_1_0707 import *
from pltMaking0707 import *
from reMatchPull0707 import *
# from stdCuts0707 import * #normal
from stdCutsMagDistCorr_out import * #trying to implement the magDistortion Correction

from getFLCdrcRefStars0707 import *
from linFLC2drc0707 import *
# from getRefiter0707 import *
# from linFLC2drcIter0707 import *
from whichIter0707 import *
# from getMatchedFLCdrc0707 import *
from getMatchedFLCdrcMagCorr2407 import *
from getZPT0707 import *
from scatterZPTdiff import *

from getRefFilt0707 import *
from linTransFilt0707 import *
from matchFiltPull0707 import *
from outEl0707 import *
from makeCMD0707 import *

# from outCorrFull import *
from magcorrTry2407 import *
    # doAll

filt_arr = ['F606W', 'F814W']
# filt_arr = ['F606W']
# filt_arr = ['F814W']

# targname_arr = ['HYDRA-II','PEGASUS-III','PHOENIX-II','RETICULUM-II','TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE','TUCANA-II-NW','TUCANA-II-SE','TUCANA-II-SW']

#,'SAGITTARIUS-II']

#,'HOROLOGIUM-I']

# targname_arr = ['BOOTES-II-NORTH','BOOTES-II-SOUTH','CARINA',
# 'CETUS-II','COLUMBA-I','DRACO-II','ERIDANUS-III','GRUS-I','GRUS-II']

# targname_arr=['HOROLOGIUM-II-EAST','HOROLOGIUM-II-WEST','INDUS-II','PICTORIS-I','PISCES-II','RETICULUM-III','SEGUE-1-EAST','SEGUE-1-WEST','SEGUE-2','SEXTANS']

# targname_arr=['TUCANA-III-EAST',
# 'TUCANA-III-WEST','TUCANA-IV-NORTH','TUCANA-IV-SOUTH','TUCANA-V','URSA-MAJOR-II-EAST','URSA-MAJOR-II-WEST','WILLMAN-1']

# targname_arr=['HOROLOGIUM-I']
# targname_arr = ['TUCANA-V']

targname_arr = ['SAGITTARIUS-II','HOROLOGIUM-I']
# targname_arr = ['PHOENIX-II']

# targname_arr = ['HYDRA-II','PEGASUS-III','PHOENIX-II','RETICULUM-II','TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE','TUCANA-II-NW','TUCANA-II-SE','TUCANA-II-SW','SAGITTARIUS-II','HOROLOGIUM-I','HOROLOGIUM-II-EAST','HOROLOGIUM-II-WEST','INDUS-II','PICTORIS-I','PISCES-II','RETICULUM-III','SEGUE-1-EAST','SEGUE-1-WEST','SEGUE-2','SEXTANS','TUCANA-III-EAST',
# 'TUCANA-III-WEST','TUCANA-IV-NORTH','TUCANA-IV-SOUTH','TUCANA-V','URSA-MAJOR-II-EAST','URSA-MAJOR-II-WEST','WILLMAN-1']


# Got all of the FLC magnitudes calculated
for c1,targname in enumerate(targname_arr):
    seDir, magCatDir, catDir = f2mag_dirs(targname,date='1305',workDir='./')
    # exp_DRC(targname)
    # outPutCorr(targname,dir=catDir)
    for c2,filt in enumerate(filt_arr):
        # get_mags(targname,filt,'1305',workDir='./')
        # wrapped(targname,filt)
        # outDiths(targname,filt,dir=catDir,suffix='_ref.dat',iter=1)
        # openFiles(targname,filt,dir=catDir,iter=1) # makes plots
        # wrapped_i(targname,filt,iter=1)
        # makeSTDcuts(catDir,targname,filt,suffix='_aftLT.dat')
        # getRef(targname,filt,dir=catDir,matchtol=50)
        # linFLC2drc(targname,filt,dir=catDir)
        # match_file = whichIter(targname,filt,dir=catDir)
        # print(match_file)
        # getMatch(targname,filt,match_file,dir=catDir,matchtol=2.5,stdTol=5)
        doAll(targname,filt,dir=catDir)
        #     doIterMatch(targname,filt,dir=catDir,matchtol=2.5,stdTol=5)
        # # doZPT(targname,filt,dir=catDir,sigTol=2.5,stdTol=0.05)
        # plotZPTdiff(targname,filt,dir=catDir,sigTol=2.5,stdTol=0.05)
        # plotMagCorrs(targname,filt,dir=catDir)

    # getRefFilt(targname,matchtol=3,dir=catDir)
    # linFiltTrans(targname,dir=catDir)
    # matchFilt(targname,dir=catDir,matchtol=3)
    # outputElen(targname, dir=catDir)
    # makeCMD(targname,dir=catDir)
