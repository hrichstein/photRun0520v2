from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import time

from hst_func import *
from f2mag import f2mag_dirs
from getJdan import getJdan
from linTrans import *

upperDir = "/Volumes/Spare Data/Hannah_Data/"
date = '1305'
# seDir, magCatDir, catDir = f2mag_dirs(date,workDir='./')
seDir, magCatDir, catDir = f2mag_dirs(date,workDir='./')

offset = 20.0 # number of pixels in one arcsecond

def distCor(targname,filt,workDir='./'):

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):
        #Load catalog
        cat = np.genfromtxt(workDir+magCatDir + jdanUse[ff]+'_'+targname+'_'+filt+'_wMag.dat',names=True)
        catCat = np.loadtxt(workDir+magCatDir + jdanUse[ff]+'_'+targname+'_'+filt+'_wMag.dat')
        #Do correction
        cor = acsDistortion(upperDir+'wfc_'+filt,cat['xr'],cat['yr'])
        #Add columns
        cat = np.hstack((catCat,cor))

        header = "flags RA DEC xr yr flux c_star magr id xc yc"
        form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f"

        np.savetxt(workDir+catDir+jdanUse[ff]+'_'+targname+'_'+filt+'_dc.dat',cat,fmt=form,header=header)

    return None

def offCor(targname,filt,workDir='./'):

    jdanUse = getJdan(targname,filt)

    if filt=='F606W':
        fils = 'f606w/'
    elif filt=='F814W':
        fils = 'f814w/'

    for jj, jdan in enumerate(jdanUse):

        # Load images, retrieve offset info
        tempim = fits.open(workDir+targname+'_'+fils+jdan+"_flc.fits")

        xoff = float(tempim[0].header["POSTARG1"])
        yoff = float(tempim[0].header["POSTARG2"])

        # Load the respective catalog
        cat = np.genfromtxt(workDir+catDir+jdan+"_"+targname+"_"+filt+"_dc.dat",names=True)
        catCat = np.loadtxt(workDir+catDir+jdan+"_"+targname+"_"+filt+"_dc.dat")

        # Create an array for the new values
        newCol = np.zeros((len(cat),2))

        # Apply offsets to columns
        newCol[:,0] = cat['xc'] - (offset * xoff)
        newCol[:,1] = cat['yc'] - (offset * yoff)

        # Combine to single array and save out
        cat = np.hstack((catCat, newCol))

        header = "flags RA DEC xr yr flux c_star magr id xc yc xo yo"
        form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f"

        np.savetxt(workDir+catDir+jdan+"_"+targname+"_"+filt+"_oc.dat",cat,header=header,fmt=form)


    return None


def matchWJCs(targname,filt,iter=0,workDir='./',matchtol=5):

    # xt, yt = 11,12
    # magr,id = 7,8

    jdanUse = getJdan(targname,filt)
    outName = "master_ids_"+targname+"_"+filt+'_{0:d}'.format(iter)+".dat"

    if iter > 0:
        suffix = "_at_{0:d}".format(iter)+".dat"
    else:
        suffix = '_oc.dat'

    master = np.genfromtxt(workDir+catDir+jdanUse[0]+"_"+targname+'_'+filt +suffix,names=True)
    masterCat = np.loadtxt(workDir+catDir+jdanUse[0]+"_"+targname+'_'+filt +suffix)

    colNs = np.array(master.dtype.names)

    if iter > 0:
        xt = np.int(np.where(colNs=='xt')[0])
        yt = np.int(np.where(colNs=='yt')[0])
        xtstr = 'xt'
        ytstr = 'yt'
    else:
        xt = np.int(np.where(colNs=='xo')[0])
        yt = np.int(np.where(colNs=='yo')[0])
        xtstr = 'xo'
        ytstr = 'yo'
    magr = np.int(np.where(colNs=='magr')[0])
    id = np.int(np.where(colNs=='id')[0])
    # Create an array of zeros with columns equal to the number of non-master dithers to store the matching id for each source
    matchids = np.zeros((len(master), (len(jdanUse)-1)))
    master = np.hstack((masterCat, matchids))

    # Loop through other images
    for dd in range(len(jdanUse)-1):
        # Load catalogs
        cat = np.genfromtxt(workDir+catDir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix,names=True)
        catCat = np.loadtxt(workDir+catDir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix)

        colCNs = np.array(cat.dtype.names)

        nF = True
        row = 0

        while (nF): # not finished
            matchrows = cat[(abs(master[row][xt] - cat[xtstr]) <= matchtol) & (abs(master[row][yt] - cat[ytstr]) <= matchtol)]

    #         # Setting the proper column number to the matching index.
            if (len(matchrows) == 1):
              master[row][xt+dd+2] = matchrows[0][id]
              row = row + 1

            elif (len(matchrows) > 1):
                magDif = np.zeros((len(matchrows),1))
                for mm in range(len(matchrows)):
                    magDif[mm] = master[row][magr] - matchrows[mm][magr]
                    small = np.argmin(magDif)
                    master[row][xt+dd+2] = matchrows[small][id]
                row += 1

            else:
              master = np.delete(master, row, 0)

            if (row >= len(master)):
                nF = False

    header =  "flags RA DEC xr yr flux c_star magr id xc yc xt yt id2 id3 id4"
    form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f %d %d %d"

    np.savetxt(workDir+catDir+outName,master, header=header, fmt=form)


    return None

# def wcsTrans(targname,filt,iter,workDir='./'):
#
#     jdanUse = getJdan(targname,filt)
#
#     if filt=='F606W':
#         fils = 'f606w/'
#     elif filt=='F814W':
#         fils = 'f814w/'
#
#     cat = np.genfromtxt(workDir+catDir+"master_ids_"+targname+"_"+filt+'_{0:d}'.format(iter)+".dat", names=True)
#     catCat = np.loadtxt(workDir+catDir+"master_ids_"+targname+"_"+filt+'_{0:d}'.format(iter)+".dat")
#
#     colNs = np.array(cat.dtype.names)
#
#     xc = np.int(np.where(colNs=='xc')[0])
#     yc = np.int(np.where(colNs=='yc')[0])
#
#     newCols = np.zeros((len(cat),2))
#
#     image = fits.open(workDir+targname+'_'+fils+jdanUse[0]+"_flc.fits")
#     w = wcs.WCS(header=image[1].header,fobj=image)
#
#     # Applied to distortion corrected values
#     # Since drawing from file, the "master" is just the first one listed.
#     newCols[:,0], newCols[:,1] = w.wcs_pix2world(catCat[:,xc],catCat[:,yc],1)
#
#     image.close()
#
#     cat = np.hstack((catCat, newCols))
#
#     header = "flags RA DEC xr yr flux c_star magr id xc yc xt yt id2 id3 id4 wcsRA wcsDEC"
#
#     form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f %d %d %d %1.7f %1.7f"
#
#
#     np.savetxt(workDir+catDir+targname+'_'+filt+'_coords_{0:d}.dat'.format(iter),cat,header=header, fmt=form)
#
#     return None

#
def pullMags(targname,filt,iter,workDir='./'):

    jdanUse = getJdan(targname,filt)

    master = np.genfromtxt(workDir+catDir+targname+'_'+filt+'_coords_{0:d}.dat'.format(iter),names=True)
    masterCat = np.loadtxt(workDir+catDir+targname+'_'+filt+'_coords_{0:d}.dat'.format(iter))

    colNs = np.array(master.dtype.names)

    if iter > 0:
        suffix = "_at_{0:d}".format(iter)+".dat"
        xt = np.int(np.where(colNs=='xo')[0])
        yt = np.int(np.where(colNs=='yo')[0])
    else:
        suffix = '_oc.dat'
        xt = np.int(np.where(colNs=='xt')[0])
        yt = np.int(np.where(colNs=='yt')[0])

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

    wra_id = np.int(np.where(colNs=='wcsRA')[0])
    wdec_id = np.int(np.where(colNs=='wcsDEC')[0])
    id = np.int(np.where(colNs=='id')[0])

    coordRows = masterCat[:,[wra_id,wdec_id,flu_id,fla_id,cs_id]]

    nCo = len(jdanUse)*int(9)
    newCols = np.zeros((len(coordRows), nCo))

    jj = 0
    cc = 0
    while jj < len(jdanUse):
        cat = np.genfromtxt(workDir+catDir+jdanUse[jj]+"_"+targname+"_"+filt+suffix,names=True)
        catCat = np.loadtxt(workDir+catDir+jdanUse[jj]+"_"+targname+"_"+filt+suffix)

        if jj==0:
            idcol = id
        else:
            idcol = xt+jj+1

        rowsMast = np.transpose(masterCat)

        newIDcol = rowsMast[idcol]
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

    header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1  xr2 yr2 xc2 yc2 xt2 yt2 xr3 yr3 xc3 yc3 xt3 yt3 xr4 yr4 xc4 xc4 xt4 yt4'
    form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f'

    np.savetxt(workDir+catDir+targname+'_'+filt+'_magList_{0:d}.dat'.format(iter),magList,header=header, fmt=form)


    return None


wcsRA, wcsDEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xc1, yc1, xt1, yt1, xr2, yr2, xc2, yc2, xt2, yt2, xr3, yr3, xc3, yc3, xt3, yt3, xr4, yr4, xc4, yc4, xt4, yt4 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40


def distStd(data,idx,idx_arr=None):

    dd = int(idx)

    if idx_arr=='arr1':
        idx_x = np.array([xt1,xt2,xt3],dtype=int)
        idx_y = np.array([yt1,yt2,yt3],dtype=int)
        x_std = np.nanstd(np.array([data[dd,idx_x[0]],data[dd,idx_x[1]],data[dd,idx_x[2]]]))
        y_std = np.nanstd(np.array([data[dd,idx_y[0]],data[dd,idx_y[1]],data[dd,idx_y[2]]]))

    elif idx_arr=='arr2':
        idx_x = np.array([xt4,xt2,xt3],dtype=int)
        idx_y = np.array([yt4,yt2,yt3],dtype=int)
        x_std = np.nanstd(np.array([data[dd,idx_x[0]],data[dd,idx_x[1]],data[dd,idx_x[2]]]))
        y_std = np.nanstd(np.array([data[dd,idx_y[0]],data[dd,idx_y[1]],data[dd,idx_y[2]]]))

    elif idx_arr=='arr3':
        idx_x = np.array([xt1,xt2,xt4],dtype=int)
        idx_y = np.array([yt1,yt2,yt4],dtype=int)
        x_std = np.nanstd(np.array([data[dd,idx_x[0]],data[dd,idx_x[1]],data[dd,idx_x[2]]]))
        y_std = np.nanstd(np.array([data[dd,idx_y[0]],data[dd,idx_y[1]],data[dd,idx_y[2]]]))

    elif idx_arr=='arr4':
        idx_x = np.array([xt1,xt3,xt4],dtype=int)
        idx_y = np.array([yt1,yt3,yt4],dtype=int)
        x_std = np.nanstd(np.array([data[dd,idx_x[0]],data[dd,idx_x[1]],data[dd,idx_x[2]]]))
        y_std = np.nanstd(np.array([data[dd,idx_y[0]],data[dd,idx_y[1]],data[dd,idx_y[2]]]))

    else:
        x_std = np.nanstd(np.array([data[dd,xt1],data[dd,xt2],data[dd,xt3],data[dd,xt4]]))
        y_std = np.nanstd(np.array([data[dd,yt1],data[dd,yt2],data[dd,yt3],data[dd,yt4]]))

    dist_std = np.sqrt(x_std**2 + y_std**2)


    return dist_std


def coordStd(data,idx,idx_arr=None):

    dd = int(idx)

    if idx_arr=='arr1':
        idx_ra = np.array([ra1,ra2,ra3],dtype=int)
        idx_dec = np.array([dec1,dec2,dec3],dtype=int)
        ra_std = np.nanstd(np.array([data[dd,idx_ra[0]],data[dd,idx_ra[1]],data[dd,idx_ra[2]]]))
        dec_std = np.nanstd(np.array([data[dd,idx_dec[0]],data[dd,idx_dec[1]],data[dd,idx_dec[2]]]))

    elif idx_arr=='arr2':
        idx_ra = np.array([ra4,ra2,ra3],dtype=int)
        idx_dec = np.array([dec4,dec2,dec3],dtype=int)
        ra_std = np.nanstd(np.array([data[dd,idx_ra[0]],data[dd,idx_ra[1]],data[dd,idx_ra[2]]]))
        dec_std = np.nanstd(np.array([data[dd,idx_dec[0]],data[dd,idx_dec[1]],data[dd,idx_dec[2]]]))

    elif idx_arr=='arr3':
        idx_ra = np.array([ra1,ra2,ra4],dtype=int)
        idx_dec = np.array([dec1,dec2,dec4],dtype=int)
        ra_std = np.nanstd(np.array([data[dd,idx_ra[0]],data[dd,idx_ra[1]],data[dd,idx_ra[2]]]))
        dec_std = np.nanstd(np.array([data[dd,idx_dec[0]],data[dd,idx_dec[1]],data[dd,idx_dec[2]]]))

    elif idx_arr=='arr4':
        idx_ra = np.array([ra1,ra3,ra4],dtype=int)
        idx_dec = np.array([dec1,dec3,dec4],dtype=int)
        ra_std = np.nanstd(np.array([data[dd,idx_ra[0]],data[dd,idx_ra[1]],data[dd,idx_ra[2]]]))
        dec_std = np.nanstd(np.array([data[dd,idx_dec[0]],data[dd,idx_dec[1]],data[dd,idx_dec[2]]]))

    else:
        ra_std = np.nanstd(np.array([data[dd,ra1],data[dd,ra2],data[dd,ra3],data[dd,ra4]]))
        dec_std = np.nanstd(np.array([data[dd,dec1],data[dd,dec2],data[dd,dec3],data[dd,dec4]]))

    coord_std = np.sqrt(ra_std**2 + dec_std**2)


    return coord_std


def filterMags(targname,filt,iter,stdCut=2.5,workDir='./'):

    data = np.genfromtxt(workDir+catDir+targname+'_'+filt+'_magList_{0:d}.dat'.format(iter))

    newCol = np.zeros((len(data),7))

    for dd in range(len(data)):
        idx = np.array([mag1,mag2,mag3,mag4],dtype=int) # columns of mags

        tmp_std_arr = np.zeros((4))

        arr1 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[2]]]) #3 is out
        arr2 = np.array([data[dd,idx[3]], data[dd,idx[1]], data[dd,idx[2]]]) #0 is out
        arr3 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[3]]]) #2 is out
        arr4 = np.array([data[dd,idx[0]], data[dd,idx[2]], data[dd,idx[3]]]) #1 is out

        mean1 = np.nanmean(arr1)
        mean2 = np.nanmean(arr2)
        mean3 = np.nanmean(arr3)
        mean4 = np.nanmean(arr4)

        std1 = np.nanstd(arr1)
        std2 = np.nanstd(arr2)
        std3 = np.nanstd(arr3)
        std4 = np.nanstd(arr4)

        tmp_std_arr[0] = abs(data[dd,idx[3]] - mean1)/std1
        tmp_std_arr[1] = abs(data[dd,idx[0]] - mean2)/std2
        tmp_std_arr[2] = abs(data[dd,idx[2]] - mean3)/std3
        tmp_std_arr[3] = abs(data[dd,idx[1]] - mean4)/std4

        max_std = np.nanmax(tmp_std_arr)

        num_above_std = (tmp_std_arr >= stdCut).sum()
        num_above_std = int(num_above_std)

        if max_std >= stdCut:
            tmp_idx = np.argsort(tmp_std_arr)[::-1][0]
            if tmp_idx == 0:
                newCol[dd,0] = mean1
                newCol[dd,1] = std1
                newCol[dd,2] = distStd(data,dd,'arr1')
                newCol[dd,3] = coordStd(data,dd,'arr1')
                newCol[dd,4] = 1 # flag cut
                newCol[dd,5] = 3 # idx cut

            elif tmp_idx == 1:
                newCol[dd,0] = mean2
                newCol[dd,1] = std2
                newCol[dd,2] = distStd(data,dd,'arr2')
                newCol[dd,3] = coordStd(data,dd,'arr2')
                newCol[dd,4] = 1
                newCol[dd,5] = 0

            elif tmp_idx == 2:
                newCol[dd,0] = mean3
                newCol[dd,1] = std3
                newCol[dd,2] = distStd(data,dd,'arr3')
                newCol[dd,3] = coordStd(data,dd,'arr3')
                newCol[dd,4] = 1
                newCol[dd,5] = 2

            else:
                newCol[dd,0] = mean4
                newCol[dd,1] = std4
                newCol[dd,2] = distStd(data,dd,'arr4')
                newCol[dd,3] = coordStd(data,dd,'arr4')
                newCol[dd,4] = 1
                newCol[dd,5] = 1

        else:
            newCol[dd,0] = np.nanmean(data[dd][idx])
            newCol[dd,1] = np.nanstd(data[dd][idx])
            newCol[dd,2] = distStd(data,dd)
            newCol[dd,3] = coordStd(data,dd)
            newCol[dd,4] = 0
            newCol[dd,5] = 4 # means no index was cut

        newCol[dd,6] = num_above_std

    data = np.hstack((data, newCol))

    header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1 xr2 yr2 xc2 yc2 xt2 yt2 xr3 yr3 xc3 yc3 xt3 yt3 xr4 yr4 xc4 xc4 xt4 yt4 mean stdev pix_std coord_std cut_flag idx_cut num_abv_std'

    form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %d %d %d'

    np.savetxt(workDir+catDir+targname+'_'+filt+'_cut_std_{0:d}.dat'.format(iter),data,header=header, fmt=form)

    return None


def addTranscols(targname,filt,iter=1,workDir='./'):

    jdanUse = getJdan(targname,filt)

    for ff in range(len(jdanUse)):

        cat = np.loadtxt(workDir+catDir+jdanUse[ff]+"_"+targname+"_"+filt+"_dc.dat")

        transCat = np.loadtxt(workDir+catDir+jdanUse[ff]+"_"+targname+"_"+filt+"_t_"+str(iter)+".dat".format(iter))

        newCol = np.zeros((len(cat),2))

        newCol[:,0] = transCat[:,0]
        newCol[:,1] = transCat[:,1]

        cat = np.hstack((cat, newCol))

        header = 'flags RA DEC xr yr flux c_star magr id xc yc xt yt'

        np.savetxt(workDir+catDir+jdanUse[ff] + "_"+targname+"_"+filt+"_at_"+str(iter)+".dat", cat, fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.6f %1.6f", header=header)


    return None


def wrapAll(targname,filt,date,workDir='./',matchtol=5,iter=1,stdCut=2.5):
#
    if iter==0:
        distCor(targname,filt,workDir=workDir)
        offCor(targname,filt,workDir=workDir)
    else:
        linTrans(targname,filt,iter,catDir)
        addTranscols(targname,filt,iter=iter,workDir='./')
#     matchWJCs(targname,filt,iter,workDir=workDir,matchtol=matchtol)
#     wcsTrans(targname,filt,iter,workDir=workDir)
#     pullMags(targname,filt,iter,workDir)
#     filterMags(targname,filt,iter,workDir=workDir,stdCut=stdCut)

    return None
