import numpy as np
import os
from astropy.io import fits
import time

from getJdan import getJdan
from f2mag0707 import f2mag_dirs
from hst_func import *
from linTrans import *

upperDir = "/Volumes/Spare Data/Hannah_Data/"
offset = 20.0

def matchWJCs_i(targname,filt,workDir='./',matchtol=0.5,iter=1):

    # xt, yt = 11,12
    # magr,id = 7,8

    jdanUse = getJdan(targname,filt)
    outName = "master_ids_"+targname+"_"+filt+"_aftLT.dat"

    master = np.genfromtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt+"_oc.dat",names=True)
    masterCat = np.loadtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt+'_oc.dat')

    colNs = np.array(master.dtype.names)

    xt = np.int(np.where(colNs=='xo')[0])
    yt = np.int(np.where(colNs=='yo')[0])
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
        cat = np.genfromtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+"_t{0:d}.dat".format(iter),names=True)
        catCat = np.genfromtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+"_t{0:d}.dat".format(iter))

        nF = True
        row = 0

        while (nF): # not finished
            matchrows = cat[(abs(master[row][xt] - cat[xtstr]) <= matchtol) & (abs(master[row][yt] - cat[ytstr]) <= matchtol)]

    #         # Setting the proper column number to the matching index.
            if (len(matchrows) == 1):
              matchids[row][dd] = matchrows[0][id]
              row = row + 1

            elif (len(matchrows) > 1):
                distDiff = np.zeros((len(matchrows),1))

                for mm in range(len(matchrows)):
                    distDiff[mm] = np.sqrt( (master[row][xt] - matchrows[mm][xtstr])**2 +  (master[row][yt] - matchrows[mm][ytstr])**2)
                    small = np.argmin(distDiff)
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

    print(targname,filt,len(master))
    np.savetxt(workDir+outName,outArr, header=header, fmt=form)


    return None


def pullMags_i(targname,filt,dir='./',suffix='_aftLT.dat',iter=1):

    jdanUse = getJdan(targname,filt)

    master = np.genfromtxt(dir+'master_ids_'+targname+'_'+filt+suffix,names=True)
    masterCat = np.loadtxt(dir+'master_ids_'+targname+'_'+filt+suffix)

    colNs = np.array(master.dtype.names)

    ra_id = np.int(np.where(colNs=='RA')[0])
    dec_id = np.int(np.where(colNs=='DEC')[0])
    flu_id = np.int(np.where(colNs=='flux')[0])
    fla_id = np.int(np.where(colNs=='flags')[0])
    cs_id = np.int(np.where(colNs=='c_star')[0])
    magr = np.int(np.where(colNs=='magr')[0])

    xr = np.int(np.where(colNs=='xr')[0])
    yr = np.int(np.where(colNs=='yr')[0])
    xc = np.int(np.where(colNs=='xc')[0])
    yc = np.int(np.where(colNs=='yc')[0])

    id2 = np.int(np.where(colNs=='id2')[0])
    id3 = np.int(np.where(colNs=='id3')[0])
    id4 = np.int(np.where(colNs=='id4')[0])

    id = np.int(np.where(colNs=='id')[0])
    coordRows = masterCat[:,[ra_id,dec_id,flu_id,fla_id,cs_id]]

    nCo = len(jdanUse)*int(9) # 4 is number of dithers
    newCols = np.zeros((len(coordRows), nCo))

    # rowsMast = np.transpose(masterCat)

    jj = 0
    cc = 0
    while jj < len(jdanUse):
        xt = np.int(np.where(colNs=='xt')[0])
        yt = np.int(np.where(colNs=='yt')[0])

        if jj==0:
            suff = '_oc.dat'
        else:
            suff = "_t{0:d}.dat".format(iter)

        cat = np.genfromtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suff,names=True)
        catCat = np.loadtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suff)

        # print('Length Master:',len(master))
        # print('Length Dither {0:d}:'.format(jj),len(cat))

        if jj==0:
            idcol = id
        elif jj==1:
            idcol = id2
            xt += 2 # the non-master dithers have two extra cols
            yt += 2
        elif jj==2:
            idcol = id3
            xt += 2
            yt += 2
        elif jj==3:
            idcol = id4
            xt += 2
            yt += 2

        newIDcol = masterCat[:,idcol]
        idx = np.asarray(newIDcol,int)

        reg = catCat[idx]

        newCols[:,cc] = reg[:,magr]
        newCols[:,cc+jj+4] = reg[:,ra_id]
        newCols[:,cc+jj+5] = reg[:,dec_id]

        newCols[:,cc+jj+12] = reg[:,xr]
        newCols[:,cc+jj+13] = reg[:,yr]

        newCols[:,cc+jj+20] = reg[:,xc]
        newCols[:,cc+jj+21] = reg[:,yc]

        newCols[:,cc+jj+28] = reg[:,xt]
        newCols[:,cc+jj+29] = reg[:,yt]

        cc += 1
        jj += 1


    magList = np.hstack((coordRows, newCols))

    header = 'RA DEC flux flags c_star mag1 mag2 mag3 mag4 ra1 '
    header += 'dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xr2 yr2 xr3 yr3 xr4 yr4 '
    header += 'xc1 yc1 xc2 yc2 xc3 yc3 xc4 yc4 xt1 yt1 xt2 yt2 xt3 yt3 xt4 yt4'

    form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f '
    form +='%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f '
    form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
    form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
    form +='%1.4f %1.4f'

    np.savetxt(dir+'matched_w_MagsPos_'+filt+suffix,magList,header=header,fmt=form)

    return None


def wrapped_i(targname,filt,iter=1):
    seDir, magCatDir, catDir = f2mag_dirs(targname,date='1305',workDir='./')
    matchWJCs_i(targname,filt,workDir=catDir,matchtol=3,iter=iter)
    pullMags_i(targname,filt,dir=catDir,suffix='_aftLT.dat',iter=iter)

    return None

    #
