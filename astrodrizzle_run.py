# Must be run in activated astroconda environment
# do conda activate astroconda, then open ipython
import numpy as np
import os
import shutil
from drizzlepac import tweakreg
from drizzlepac import astrodrizzle
from makeFLCdirs import run_command

from drizzlepac import staticMask, sky, adrizzle, createMedian, ablot, drizCR

upperDir = '/Volumes/Spare Data/Hannah_Data/'

def run_aDriz(targnames):

    names = np.loadtxt(targnames,dtype=str)

    for nn in range(len(names)):
        temp_dir_606 = names[nn] + '_f606w/'
        temp_dir_814 = names[nn] + '_f814w/'

        astrodrizzle.AstroDrizzle(temp_dir_606 + "*.fits",\
                          configobj=upperDir+'astrodrizzle_new.cfg')
        #
        astrodrizzle.AstroDrizzle(temp_dir_814 + "*.fits",\
                          configobj=upperDir+'astrodrizzle_new.cfg')

    return None
