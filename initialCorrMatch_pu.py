import numpy as np
import os
from astropy.io import fits

# from getJdan import getJdan
# from f2mag0707 import f2mag_dirs
from hst_func import *
from linTrans import *

upperDir = "/Volumes/Spare Data/Hannah_Data/"
offset = 20.0 # number of pixels in one arcsecond
# seDir, magCatDir, catDir = f2mag_dirs(targname,date='1305',workDir='./')

# Eventually put the magCorrection function in here

def distCor(targname,filt,jdanUse,workDir='./'):

    # jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):
        #Load catalog
        cat = np.genfromtxt(workDir+ jdanUse[ff]+'_'+filt+'photU.dat',names=True)
        catCat = np.loadtxt(workDir+ jdanUse[ff]+'_'+filt+'photU.dat')

        colNs = np.array(cat.dtype.names)
        #Do correction
        cor = acsDistortion(upperDir+'wfc_'+filt,cat['xcenter'],cat['ycenter'])
        #Add columns
        cat = np.hstack((catCat,cor))

        s0 = ' '
        header = s0.join(colNs)
        header += ' x_dc y_dc'

        np.savetxt(workDir+jdanUse[ff]+'_'+targname+'_'+filt+'_dc.dat',cat,header=header)


    return None


def offCor(targname,filt,jdanUse,workDir='./'):

    # jdanUse = getJdan(targname,filt)

    if filt=='F606W':
        fils = 'f606w/'
    elif filt=='F814W':
        fils = 'f814w/'

    for jj, jdan in enumerate(jdanUse):

        # Load images, retrieve offset info
        tempim = fits.open(targname+'_'+fils+jdan+"_flc.fits")

        xoff = float(tempim[0].header["POSTARG1"])
        yoff = float(tempim[0].header["POSTARG2"])

        tempim.close()

        # Load the respective catalog
        cat = np.genfromtxt(workDir+jdan+"_"+targname+"_"+filt+"_dc.dat",names=True)
        catCat = np.loadtxt(workDir+jdan+"_"+targname+"_"+filt+"_dc.dat")

        colNs = np.array(cat.dtype.names)

        # Create an array for the new values
        newCol = np.zeros((len(cat),2))

        # Apply offsets to columns
        newCol[:,0] = cat['x_dc'] - (offset * xoff)
        newCol[:,1] = cat['y_dc'] - (offset * yoff)

        # Combine to single array and save out
        cat = np.hstack((catCat, newCol))

        s0 = ' '
        header = s0.join(colNs)
        header += ' x_oc y_oc'

        np.savetxt(workDir+jdan+"_"+targname+"_"+filt+"_oc.dat",cat,header=header)


    return None

# Specifically for first run through
def matchWJCs(targname,filt,jdanUse,workDir='./',matchtol=2,suffix='_ref.dat'):

    # xt, yt = 11,12
    # magr,id = 7,8

    # jdanUse = getJdan(targname,filt)
    outName = "master_ids_"+targname+"_"+filt+suffix

    master = np.genfromtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt+"_oc.dat",names=True)
    masterCat = np.loadtxt(workDir+jdanUse[0]+"_"+targname+'_'+filt+'_oc.dat')

    colNs = np.array(master.dtype.names)

    xt = np.int(np.where(colNs=='x_oc')[0])
    yt = np.int(np.where(colNs=='y_oc')[0])
    xtstr = 'x_oc'
    ytstr = 'y_oc'

    magr = np.int(np.where(colNs=='magr')[0])
    id = np.int(np.where(colNs=='id')[0])
    # Create an array of zeros with columns equal to the number of non-master dithers to store the matching id for each source

    # master = np.hstack((masterCat, matchids))

    # Loop through other images
    for dd in range(len(jdanUse)-1):
        # Load catalogs
        cat = np.genfromtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+'_oc.dat',names=True)
        catCat = np.loadtxt(workDir+jdanUse[dd+1]+"_"+targname+'_'+filt+'_oc.dat')

        # print('Length Cat',len(cat))

        # colNs = np.array(cat.dtype.names)

        nF = True
        row = 0

        while (nF): # not finished
            matchrows = cat[(abs(master[row][xtstr] - cat[xtstr]) <= matchtol) & (abs(master[row][ytstr] - cat[ytstr]) <= matchtol)]

    # Setting the proper column number to the matching index.
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
                # print('Max MatchID',np.max(matchids))

    outArr = np.hstack((masterCat,matchids))

    s0 = ' '
    header = s0.join(colNs)
    header +=  " id2 id3 id4"

    print(targname,filt,len(master))
    np.savetxt(workDir+outName,outArr,header=header)


    return None


def pullMags(targname,filt,jdanUse,dir='./',suffix='_0.dat'):

    # jdanUse = getJdan(targname,filt)

    master = np.genfromtxt(dir+'master_ids_'+targname+'_'+filt+suffix,names=True)
    masterCat = np.loadtxt(dir+'master_ids_'+targname+'_'+filt+suffix)

    colNs = np.array(master.dtype.names)

    flu_id = np.int(np.where(colNs=='final_phot')[0])
    magr = np.int(np.where(colNs=='magr')[0])

    xr = np.int(np.where(colNs=='xcenter')[0])
    yr = np.int(np.where(colNs=='ycenter')[0])
    xc = np.int(np.where(colNs=='x_dc')[0])
    yc = np.int(np.where(colNs=='y_dc')[0])

    id2 = np.int(np.where(colNs=='id2')[0])
    id3 = np.int(np.where(colNs=='id3')[0])
    id4 = np.int(np.where(colNs=='id4')[0])
    # xt, yt are xo,yo in dithers 2-4 as the transformed positions
    # don't change for the first dither
    xt = np.int(np.where(colNs=='x_oc')[0])
    yt = np.int(np.where(colNs=='y_oc')[0])

    id = np.int(np.where(colNs=='id')[0])
    coordRows = masterCat[:,[flu_id]]

    nCo = len(jdanUse)*int(7) # 4 is number of dithers
        # one col for mag, 6 for coordinates
    newCols = np.zeros((len(coordRows), nCo))

    # rowsMast = np.transpose(masterCat)

    jj = 0
    cc = 0
    while jj < len(jdanUse):
        suff = '_oc.dat'
        cat = np.genfromtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suff,names=True)
        catCat = np.loadtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suff)

        if jj==0:
            idcol = id
        elif jj==1:
            idcol = id2
        elif jj==2:
            idcol = id3
        elif jj==3:
            idcol = id4


        newIDcol = masterCat[:,idcol]
        idx = np.asarray(newIDcol,int)

        # print('Cat Length',len(cat))
        # print('Max IDX',np.max(idx))
    #
        reg = catCat[idx]

        newCols[:,cc] = reg[:,magr]

        newCols[:,cc+jj+4] = reg[:,xr]
        newCols[:,cc+jj+5] = reg[:,yr]

        newCols[:,cc+jj+12] = reg[:,xc]
        newCols[:,cc+jj+13] = reg[:,yc]

        newCols[:,cc+jj+20] = reg[:,xt]
        newCols[:,cc+jj+21] = reg[:,yt]
    #
        cc += 1
        jj += 1
    #
    #
    magList = np.hstack((coordRows, newCols))

    header = 'flux mag1 mag2 mag3 mag4 xr1 yr1 xr2 yr2 xr3 yr3 xr4 yr4 xc1 yc1 xc2 yc2 xc3 yc3 xc4 yc4 xt1 yt1 xt2 yt2 xt3 yt3 xt4 yt4'

    np.savetxt(dir+'matched_w_MagsPos'+suffix,magList,header=header)

    return None


def wrapped(targname,filt,jdan=None,catDir='./'):
    distCor(targname,filt,jdan,workDir=catDir)
    offCor(targname,filt,jdan,workDir=catDir)
    matchWJCs(targname,filt,jdan,workDir=catDir,matchtol=3,suffix='_0.dat')
    pullMags(targname,filt,jdan,dir=catDir,suffix='_0.dat')
    matchWJCs(targname,filt,jdan,workDir=catDir,matchtol=0.5,suffix='_ref.dat')
    pullMags(targname,filt,jdan,dir=catDir,suffix='_ref.dat')

    return None
    #
