import numpy as np
import os
from astropy.io import fits

from getDRCfiltRef_pu import *
    # getRefDRCFilt
# from drcFiltLinTrans_pu import *
    # linFiltTransDRC
# from matchDRCfilt_pu import *
    # matchFiltDRC

# from getJdan import getJdan
from runPU_1008 import f2mag_dirs, runPhotUtils
from initialCorrMatch_pu import *
    # distCor, offCor, matchWJCs, pullMags, wrapped
# from linTrans_1_pu import outDiths, makePlot, openFiles
from reMatchPull_pu import *
    # matchWJCs_i, pullMags_i, wrapped_i
# from stdCuts_pu import makeSTDcuts
# from drcFLCref_pu import *
    # getRef
# from linFLC2drc_pu import *
    # linFLC2drc
# from whichIter_pu import whichIter, doIterMatch, getMatch
# from getZPT_pu import *
#     # getZPT, applyZPT, doZPT
# from plotZPTdiff_pu import plotZPTdiff
# from getRefFilt_pu import *
#     # getRefFilt
# from linTransFiltFLC_pu import *
# from matchFiltFLC_pu import *
#     # matchFilt
# from makeCMD_pu import makeCMD
# from matchFLCdrcAll import *
#     # doIterMatchDRC
#
from make9plots import *
#     # feedFunc
# from matchOldNew import *
    # getONRef, getONmatch
# from linONtrans import *
# from makeCMDon import *
# from getRefPSFflc_pu import *
# from linPSFtoFLC_1408 import *
# from matchPSFflc_1408 import *
# from getRefAPERflc_1408 import *
# from linAPER2flc_1408 import *
# from matchAPERflc_1408 import *
# from makeCMDaperF_1408 import *
# from makeCMDpsfF_1408 import *
# from drcFLC_diff_1408bins import *
# from getRefDRCon import *
# from linDRCon_1408 import *
# from matchONdrc_1408 import *
# from matchONnoZPT_flc import *
from masterFeed_matched import *

filt_arr = ['F606W','F814W']
# filt_arr = ['F606W']
# targname_arr = ['SAGITTARIUS-II','HOROLOGIUM-I']

targname_arr = ['HOROLOGIUM-I']

for c1,targname in enumerate(targname_arr):
    rand_tuple = f2mag_dirs(targname,date='10Aug',workDir='./')
    saveDir = rand_tuple[-1]
    print(saveDir)

    # getRefDRCFilt(targname,dir=saveDir,matchtol=3)
    # linFiltTransDRC(targname,dir=saveDir)
    # matchFiltDRC(targname,dir=saveDir,matchtol=3)

    # for c2,filt in enumerate(filt_arr):
    #     jdan = getJdan(targname,filt)
        # runPhotUtils(targname,filt,jdan,saveDir=saveDir)
        # wrapped(targname,filt,jdan,catDir=saveDir)
        # outDiths(targname,filt,jdan,dir=saveDir,suffix='_ref.dat',iter=1)
        # openFiles(targname,filt,jdan,dir=saveDir,iter=1)
        # wrapped_i(targname,filt,jdan,iter=1,catDir=saveDir)
        # makeSTDcuts(saveDir,filt,suffix='_aftLT.dat')
        # getRef(targname,filt,dir=saveDir,matchtol=50)
        # linFLC2drc(targname,filt,dir=saveDir)
        ## match_file = whichIter(targname,filt,dir=saveDir)
        ## print(match_file)
        ## getMatch(targname,filt,match_file,dir=saveDir,matchtol=2.5,stdTol=5)
        # doIterMatch(targname,filt,dir=saveDir,matchtol=2.5,stdTol=5)
        # doZPT(targname,filt,dir=saveDir,sigTol=2.5,stdTol=0.05)
        # plotZPTdiff(targname,filt,dir=saveDir,sigTol=2.5,stdTol=0.05)
        # plotZPTdiff(targname,filt,dir=saveDir,sigTol=2.5,stdTol=0.05)
    # getRefFilt(targname,matchtol=3,dir=saveDir)
    # linFiltTrans(targname,dir=saveDir)
    # matchFilt(targname,dir=saveDir,matchtol=3)
    # makeCMD(targname,dir=saveDir)
    # doIterMatchDRC(targname,filt='F606W',dir=saveDir,matchtol=3,stdTol=5)
    # doIterMatchDRC(targname,filt='F814W',dir=saveDir,matchtol=3,stdTol=5)
    # feedFunc(targname,dir=saveDir)
    # getONRef(targname,filt='F606W',dir=saveDir,matchtol=10)
    # linONtrans(targname,filt='F606W',dir=saveDir)
    # getONmatch(targname,filt='F606W',dir=saveDir,matchtol=3,stdTol=5)
    # makeCMD(targname,dir=saveDir,filt='F606W')
    # getRefPSFflc('HOROLOGIUM-I',matchtol=1,dir=saveDir)
    # linPSFTrans(targname,dir=saveDir)
    # matchPSFflc(targname,dir=saveDir,matchtol=3)
    # feedFunc(targname,dir=saveDir,name='psfFLC')
    # getRefaperflc(targname,matchtol=3,dir=saveDir)
    # linAPERTrans(targname,dir=saveDir)
    # matchAPERflc(targname,dir=saveDir,matchtol=3)
    # feedFunc(targname,dir=saveDir,name='aperFLC')
    # makeCMD(targname,dir=saveDir)
    # getRefdrcON(targname,matchtol=3,dir=saveDir)
    # linDRConTrans(targname,dir=saveDir)
    # matchONdrc(targname,dir=saveDir,matchtol=3)
    # feedFunc(targname,dir=saveDir,name='flcDRCnoZPT')
    masterFeed(targname,dir=saveDir)
