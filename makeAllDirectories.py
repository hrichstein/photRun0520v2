import numpy as np
import os


def f2mag_dirs(targname,workDir='./catMatchFLCdrc18Oct'):

    # Gives the name of the directory with se files
    # Creates folder for raw mags if it doesn't exists
    # returns the names of both

    magCatDir = workDir + '/'  # + 'photUtils' + date
    catDir = magCatDir + 'catDir_' + targname + '/'

    # if not os.path.exists(os.path.join(".",magCatDir)):
    #     os.makedirs(magCatDir)
    if not os.path.exists(os.path.join(".",catDir)):
        os.makedirs(catDir)

    return None


targname_arr = np.genfromtxt('targnamesDirections2.txt',dtype='str')


def main():

    for c1,targname in enumerate(targname_arr):
        f2mag_dirs(targname,workDir='./catMatchFLCdrc18Oct')


if __name__ == '__main__':
    main()
