from astropy.io import fits
import numpy as np
import os

from getJdan import getJdan


def matchWJCs(targname,filt,workDir='catRawMags1305/catDir/',matchtol=3):

    # xt, yt = 11,12
    # magr,id = 7,8

    jdanUse = getJdan(targname,filt)
    outName = "master_ids_"+targname+"_"+filt+"_at_2705.dat"

    suffix = "_at_2705.dat"

    master = np.genfromtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt +suffix,names=True)
    masterCat = np.loadtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt +suffix)

    colNs = np.array(master.dtype.names)


    xt = np.int(np.where(colNs=='xt')[0])
    yt = np.int(np.where(colNs=='yt')[0])
    xtstr = 'xt'
    ytstr = 'yt'

    magr = np.int(np.where(colNs=='magr')[0])
    id = np.int(np.where(colNs=='id')[0])
    # Create an array of zeros with columns equal to the number of non-master dithers to store the matching id for each source
    matchids = np.zeros((len(master), (len(jdanUse)-1)))
    # master = np.hstack((masterCat, matchids))

    # Loop through other images
    for dd in range(len(jdanUse)-1):
        # Load catalogs
        cat = np.genfromtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix,names=True)
        catCat = np.loadtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix)

        nF = True
        row = 0

        while (nF): # not finished
            matchrows = cat[(abs(master[row][xtstr] - cat[xtstr]) <= matchtol) & (abs(master[row][ytstr] - cat[ytstr]) <= matchtol)]

    #         # Setting the proper column number to the matching index.
            if (len(matchrows) == 1):
              matchids[row][dd] = matchrows[0][id]
              row = row + 1

            elif (len(matchrows) > 1):
                magDif = np.zeros((len(matchrows),1))
                for mm in range(len(matchrows)):
                    magDif[mm] = master[row][magr] - matchrows[mm][magr]
                    small = np.argmin(magDif)
                    matchids[row][dd] = matchrows[small][id]
                row += 1

            else:
              master = np.delete(master, row, 0)
              masterCat = np.delete(masterCat, row, 0)
              matchids = np.delete(matchids,row,0)

            if (row >= len(master)):
                nF = False

    outArr = np.hstack((masterCat,matchids))

    header =  "flags RA DEC xr yr flux c_star magr id xc yc xt yt id2 id3 id4"
    form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f %d %d %d"

    print(len(master))
    np.savetxt(workDir+outName,outArr, header=header, fmt=form)


    return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'
matchWJCs(targname,filt,workDir='catRawMags1305/catDir/',matchtol=5)
