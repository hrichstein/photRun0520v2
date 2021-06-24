import numpy as np
import os
from astropy.io import fits

from getDRCfiltRef_pu import *
#     # getRefDRCFilt
from drcFiltLinTrans_pu import *
#     # linFiltTransDRC
from matchDRCfilt_pu import *
    # matchFiltDRC

from getJdan import getJdan
# from runPU_1008 import runPhotUtils
from runPU_3pix import runPhotUtils
# from runPU_drc import f2mag_dirs
from initialCorrMatch_pu import *
    # distCor, offCor, matchWJCs, pullMags, wrapped
from linTrans_1_pu import outDiths, makePlot, openFiles
from reMatchPull_pu import *
    # matchWJCs_i, pullMags_i, wrapped_i
from stdCuts_pu import makeSTDcuts
from drcFLCref_pu import *
    # getRef
from linFLC2drc_pu import *
    # linFLC2drc
from whichIter_pu import whichIter, doIterMatch, getMatch
# from getZPT_pu import *
#     # getZPT, applyZPT, doZPT
# from newZPT_pu2 import *
    # uses TopCat cut files
# from getZPT_pu3 import *
from newZPT_0510 import *
# from dummyZPT import *
# from plotZPTdiff_pu import plotZPTdiff
from getRefFilt_pu import *
#     # getRefFilt
from linTransFiltFLC_pu import *
from matchFiltFLC_pu import *
#     # matchFilt
# from makeCMD_pu import makeCMD
# from matchFLCdrcAll import *
#     # doIterMatchDRC
# from applyRedDm import applyRedDm
# from makeCMD_abs import makeCMDabs
# from make9plots import *
#     # feedFunc
# from drcFLC_diff_1408bins import *
# from match4cstar import match4cStar606, match4cStar814
# from match4cstar_f606wWmag import match4cStar606
# targname_arr = np.genfromtxt('targnamesDirections2.txt',dtype='str')
targname_arr = np.genfromtxt('targnamesPost.txt',dtype='str')

# targname_arr = ['TUCANA-II-SE']

filt_arr = ['F814W','F606W']
# filt_arr = ['F814W']
# targname_arr = ['SAGITTARIUS-II','HOROLOGIUM-I']

# targname_arr = ['RETICULUM-III']


def f2mag_dirs(targname,date='28Sep',workDir='./'):

    magCatDir = workDir + '/' + 'photUtils' + date + '/'
    catDir = magCatDir + 'catDir_' + targname + '/'

    if not os.path.exists(os.path.join(".",magCatDir)):
        os.makedirs(magCatDir)
    if not os.path.exists(os.path.join(".",catDir)):
        os.makedirs(catDir)

    # return workDir + 'catRawMags' + date + '/catDir_' + targname + '/'

    return catDir


for c1,targname in enumerate(targname_arr):
    saveDir = f2mag_dirs(targname,date='21Oct',workDir='.')
    # saveDir = saveDir[-1]
    print(targname)
    # getRefDRCFilt(targname,dir=saveDir,matchtol=3)
    # linFiltTransDRC(targname,dir=saveDir)
    # matchFiltDRC(targname,dir=saveDir,matchtol=3)
    getRefDRCFilt(targname,dir='./photUtils21Oct/catDir_'+targname+'/',matchtol=3)
    linFiltTransDRC(targname,dir='./photUtils21Oct/catDir_'+targname+'/')
    matchFiltDRC(targname,dir='./photUtils21Oct/catDir_'+targname+'/',matchtol=3)

    # for c2,filt in enumerate(filt_arr):
        # jdan = getJdan(targname,filt)
        # runPhotUtils(targname,filt,jdan,saveDir=saveDir)
        # wrapped(targname,filt,jdan,catDir=saveDir)
        # outDiths(targname,filt,jdan,dir=saveDir,suffix='_ref.dat',iter=1)
        # openFiles(targname,filt,jdan,dir=saveDir,iter=1)
        # wrapped_i(targname,filt,jdan,iter=1,catDir=saveDir)
        # makeSTDcuts(saveDir,filt,suffix='_aftLT.dat')
        # getRef(targname,filt,dir=saveDir,matchtol=50)
        # linFLC2drc(targname,filt,dir=saveDir)
        # # ### Below are included in doIterMatch
        # match_file = whichIter(targname,filt,dir=saveDir)
        # # ### print(match_file)
        # getMatch(targname,filt,match_file,dir=saveDir,matchtol=2.5,stdTol=5)
        # doIterMatch(targname,filt,dir=saveDir,matchtol=2.5,stdTol=5)
        # print(filt)
        # doZPT(targname,filt,dir=saveDir,sigTol=3.5,stdTol=0.1)
        # plotZPTdiff(targname,filt,dir=saveDir,sigTol=2.5,stdTol=0.05)
    #
    # getRefFilt(targname,matchtol=3,dir=saveDir)
    # getRefFilt(targname,matchtol=3,dir='./photUtils28Sep/catDir_'+targname+'/')
    # linFiltTrans(targname,dir='./photUtils28Sep/catDir_'+targname+'/',newdir=saveDir)
    # matchFilt(targname,dir='./photUtils28Sep/catDir_'+targname+'/',newdir=saveDir,matchtol=3)


    # getRefFilt(targname,matchtol=3,dir='./catRawMags20Aug/catDir_' + targname + '/')
    # linFiltTrans(targname,dir='./catRawMags20Aug/catDir_' + targname + '/',newdir=saveDir)
    # matchFilt(targname,dir='./catRawMags20Aug/catDir_' + targname + '/',newdir=saveDir,matchtol=3)
    # makeCMD(targname,dir=saveDir)
    # applyRedDm(targname,dir=saveDir)
    # makeCMDabs(targname,dir=saveDir)
    # match4cStar606(targname,dir=saveDir,matchtol=3)
    # match4cStar814(targname,dir=saveDir,matchtol=3)
    # doIterMatchDRC(targname,filt='F606W',dir=saveDir,matchtol=3,stdTol=5)
    # doIterMatchDRC(targname,filt='F814W',dir=saveDir,matchtol=3,stdTol=5)
    # feedFunc(targname,dir=saveDir)
    #
