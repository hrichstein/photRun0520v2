"""
Gets the identification information for the FLC files

"""
import numpy as np


def getJdan(targname,filt,dir='./',drc=False):

    jdanList_file = np.loadtxt(dir + targname + "_flcs.txt",dtype="<U50")

    jdan = jdanList_file[:,0]
    jdanFilter = jdanList_file[:,1]

    jdanUse = []
    for jj in range(len(jdan)):
        if jdanFilter[jj]==filt:
            jdanUse.append(jdan[jj])

    jdanUse = np.array(jdanUse)

    # The name of the corresponding DRC file
    drcUse = jdanUse[0][:6] + '010_drc.fits'

    if drc:
        return jdanUse, drcUse
    else:
        return jdanUse
