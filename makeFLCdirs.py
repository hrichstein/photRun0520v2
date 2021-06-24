# Creates two directories for each pointing, one for each filter
# Creates inner directories to house crclean fits
# Moves flc_fits files to new directories

import numpy as np
import os

upperDir = '/Volumes/Spare Data/Hannah_Data/'

def make_dir(targnames):

    names = np.loadtxt(targnames,dtype=str)

    f606w_dirs = [[] for xx in range(len(names))]
    f814w_dirs = [[] for xx in range(len(names))]

    for nn in range(len(names)):
        f814w_dirs[nn] = names[nn] + '_f814w'
        f606w_dirs[nn] = names[nn] + '_f606w'

    for dir606w in f606w_dirs:
        if not os.path.exists(os.path.join(".",dir606w)):
            os.makedirs(dir606w)
        if not os.path.exists(os.path.join(".",dir606w,"crClean")):
            os.makedirs(os.path.join(".",dir606w,"crClean"))

    for dir814w in f814w_dirs:
        if not os.path.exists(os.path.join(".",dir814w)):
            os.makedirs(dir814w)
        if not os.path.exists(os.path.join(".",dir814w,"crClean")):
            os.makedirs(os.path.join(".",dir814w,"crClean"))

    return None

def run_command(com):

    s0 = ''
    com = s0.join(com)
    res = os.system(com)

    return res

def move_files(targnames):

    names = np.loadtxt(targnames,dtype=str)

    # file_names = [[] for xx in range(len(names))]

    for nn in range(len(names)):
        file_names =  upperDir + names[nn] + '_flcs.txt'

        temp_dat = np.genfromtxt(file_names,dtype=str)

        temp_dir_814 = names[nn] + '_f814w'
        temp_dir_606 = names[nn] + '_f606w'

        for tt in range(len(temp_dat)):

            if temp_dat[tt,1] == 'F606W':
                jdan_name = temp_dat[tt,0] + ''

                com = ["mv flc_fits/", jdan_name, "_flc.fits ", temp_dir_606]
            elif temp_dat[tt,1] == 'F814W':
                jdan_name = temp_dat[tt,0] + ''

                com = ["mv flc_fits/", jdan_name, "_flc.fits ", temp_dir_814]

            run_command(com)

    return None

# def fortCode(targnames):
# # Try appending the paths to a list and using os.walk
# have a code entirely dedicated to executing the file
# then have a second code that is doing the searching for files and you could move or copy the executing code and import it to the first one to run the execute function

# fortran_key = ".f"
# directory = "." #stands for current directory
# fortran_files = []
#
# for root, dirs, files, in os.walk(".",topdown=False):
#   for file in files:
#     if(file.endswith(fortran_key)):
#       fortran_files.append(os.path.join(root,file))
#execute.py
# def execute(file):
  #do something that executes
  #main code.py
# import execute
# import shutil
#
# #searching thing and once you find the file at some path
# os.copy("execute.py",path)
# execute(file)

# run whatever os.system thing to be discovered later on fortran_files
#     names = np.loadtxt(targnames,dtype=str)
#
#     for nn in range(len(names)):
#         temp_dir_606 = names[nn] + '_f606w/crClean/'
#         temp_dir_814 = names[nn] + '_f814w/crClean/'
#
#         # com = ["cp flt2wj2_acswfc.e ",temp_dir_606]
#         # run_command(com)
#         com_1 = ["./",temp_dir_606,"./flt2wj2_acswfc.e ","jdan", "*.fits"]
#         run_command(com_1)
#
#         # com2 = ["cp flt2wj2_acswfc.e ",temp_dir_814]
#         # run_command(com2)
#         com2_1 = [temp_dir_814,"./flt2wj2_acswfc.e ","jdan", "*crclean.fits"]
#         run_command(com2_1)
#
#     return None






# End
