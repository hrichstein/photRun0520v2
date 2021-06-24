import numpy as np
# from getJdan import getJdan
from f2mag import get_mags
from f2mag import f2mag_dirs
from flcFns import *
from pMatch import *
from pRefStars import *
from linTrans import *

workDir = './'
upperDir = "/Volumes/Spare Data/Hannah_Data/"

filt_arr = ['F606W', 'F814W']
# filt_arr = ['F814W']
targname = 'HOROLOGIUM-I'

matchtol = 5.0
pixtol = 2
date = '1305'
catDir = workDir+'catRawMags'+date+'/catDir/'

# iter = int(2)
nF = True
iter = 1

while (nF):
    if iter == 0:
        for ff, filt in enumerate(filt_arr):
            get_mags(targname,filt,date,workDir=workDir)
            wrapAll(targname,filt,date,workDir=workDir,matchtol=matchtol,iter=iter,\
            stdCut=2.5)
        pMatch(targname,iter,catDir,matchtol=pixtol)
        for ff, filt in enumerate(filt_arr):
            pRefStars(targname,filt,(iter),catDir,magHi=24.5,magLo=21,\
            stdTol=0.1,posTol=5)

    else:
        for ff, filt in enumerate(filt_arr):
            # pRefStars(targname,filt,iter,catDir,magHi=24.5,magLo=21,\
            # stdTol=0.1,posTol=5)

            wrapAll(targname,filt,date,workDir=workDir,matchtol=matchtol,iter=iter,\
                stdCut=2.5)

        pMatch(targname,iter,catDir,matchtol=pixtol)
        # for ff, filt in enumerate(filt_arr):
        #     new = np.loadtxt(catDir+targname+'_refStars_'+filt+'_{0:d}.dat'.format(iter))
        #     old = np.loadtxt(catDir+targname+'_refStars_'+filt+'_{0:d}.dat'.format(iter-1))
        #
        #     if len(new)==len(old):
        #         nF = False
        #     else:
        #         iter += 1
        new = np.loadtxt(catDir+targname+'_refStars_'+'F814W'+'_{0:d}.dat'.format(iter))
        old = np.loadtxt(catDir+targname+'_refStars_'+'F814W'+'_{0:d}.dat'.format(iter-1))

        if len(new)==len(old):
            nF = False
        else:
            iter += 1
