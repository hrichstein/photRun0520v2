import numpy as np

# upperDir = "/Volumes/Spare Data/Hannah_Data/"
dir = './'

def getJdan(targname,filt):

    jdanList_file = np.loadtxt(dir+targname+"_flcs.txt",dtype="<U50")

    jdan = jdanList_file[:,0]
    jdanFilter = jdanList_file[:,1]

    jdanUse = []
    for jj in range(len(jdan)):
        if jdanFilter[jj]==filt:
            jdanUse.append(jdan[jj])

    jdanUse = np.array(jdanUse)

    return jdanUse
