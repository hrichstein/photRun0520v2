import shutil
import os
import numpy as np

def mv_crclean(targnames):

    names = np.loadtxt(targnames,dtype=str)

    for nn in range(len(names)):
        temp_dir_606 = names[nn] + '_f606w/'
        temp_dir_814 = names[nn] + '_f814w/'

        files = os.listdir(temp_dir_606)
        files2 = os.listdir(temp_dir_814)

        for f in files:
            if (f.endswith("_crclean.fits")):
                shutil.move(temp_dir_606+f, temp_dir_606+'crClean')

        for f in files2:
            if (f.endswith("_crclean.fits")):
                shutil.move(temp_dir_814+f, temp_dir_814+'crClean')

    return None
