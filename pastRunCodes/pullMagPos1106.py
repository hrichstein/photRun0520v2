import numpy as np
from getJdan import *

def pullMags(targname,filt,dir='./',suffix='_at_1106.dat'):

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
    # xt, yt are xo,yo in dithers 2-4 as the transformed positions
    # don't change for the first dither
    xt = np.int(np.where(colNs=='xt')[0])
    yt = np.int(np.where(colNs=='yt')[0])

    id = np.int(np.where(colNs=='id')[0])
    coordRows = masterCat[:,[ra_id,dec_id,flu_id,fla_id,cs_id]]

    nCo = len(jdanUse)*int(9) # 4 is number of dithers
    newCols = np.zeros((len(coordRows), nCo))

    # rowsMast = np.transpose(masterCat)

    jj = 0
    cc = 0
    while jj < len(jdanUse):
        suffix = '_at_2705.dat'
        cat = np.genfromtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suffix,names=True)
        catCat = np.loadtxt(dir+jdanUse[jj]+"_"+targname+"_"+filt+suffix)

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

    np.savetxt(dir+'matched_w_MagsPos1106.dat',magList,header=header,fmt=form)

    return None

targname='HOROLOGIUM-I'
filt='F606W'
dir = 'catRawMags1305/catDir/'

pullMags(targname,filt,dir=dir,suffix='_at_1106.dat')
