import numpy as np
import os
# import matplotlib.pyplot as plt

def matchFiltDRC(targname,dir='./',matchtol=3):

    x = os.listdir(dir)
    for ii in x:
        if ii.endswith('_F606WphotU.dat'):
            f606w_file = ii
        elif ii.endswith('_F814WphotU.dat'):
            f814w_file = ii

    f606wN = np.genfromtxt(dir+f606w_file,names=True)
    f606w = np.genfromtxt(dir+f606w_file)

    f814wN = np.genfromtxt(dir+
    targname+'_DRCfiltTrans2_pU.dat',names=True)
    f814w = np.genfromtxt(dir+targname+'_DRCfiltTrans2_pU.dat')

    col606 = np.array(f606wN.dtype.names)
    col814 = np.array(f814wN.dtype.names)

    # Getting columns of x,y of transformed F814W to F606W
    xI = np.int(np.where(col814=='x_f606wTrans')[0])
    yI = np.int(np.where(col814=='y_f606wTrans')[0])
    magI = np.int(np.where(col814=='magr')[0])
    idI = np.int(np.where(col814=='id')[0])

    xV = np.int(np.where(col606=='xcenter')[0])
    yV = np.int(np.where(col606=='ycenter')[0])
    magV = np.int(np.where(col606=='magr')[0])
    idV = np.int(np.where(col606=='id')[0])

    master_in = f606w[:,[idV,xV,yV,magV]]
    idV_mas,x,y,magrV = 0,1,2,3

    cat = f814w
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    len606 = len(f606w)
    len814 = len(f814w)

    minLen = np.min([len606,len814])

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xI,yI,idI)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 1
            if matchtol <= 10:
                master_in = f606w[:,[idV,xV,yV,magV]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen)

    idV_mas, xV_mas, yV_mas, magV_mas, idI_mas = 0, 1, 2, 3, 4

    # newCols = np.zeros((len(master),3))

    idCol606 = master[:,idV_mas]
    idx606 = np.asarray(idCol606,int)
    reg606 = f606w[idx606]

    idCol814 = master[:,idI_mas]
    idx814 = np.asarray(idCol814,int)
    reg814 = f814w[idx814]

    outArr = np.hstack((reg606,reg814))
    header = 'id_f606w xcenter_f606w ycenter_f606w aperture_sum_f606w annulus_median_f606w aper_bkg_f606w final_phot_f606w magr_f606w id_f814w xcenter_f814w ycenter_f814w aperture_sum_f814w annulus_median_f814w aper_bkg_f814w final_phot_f814w magr_f814w x_f606wTrans y_f606wTrans'

    np.savetxt(dir+targname+'_filtMatchDRC_pU.dat',outArr,header=header)


    return None

def matchlistID(master,cat,matchtol,x1,y1,x2,y2,id_mat):

    matchids_in = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:

        matchrows = cat[(abs(master[row][x1] - cat[:,x2]) \
            <= matchtol) & (abs(master[row][y1] - cat[:,y2])<= matchtol)]

        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[0][id_mat]
            row += 1

        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt( (master[row][x1] - \
                matchrows[dd][x2])**2 +  (master[row][y1] \
                                    - matchrows[dd][y2])**2)
            small = np.argmin(distDiff)
            matchids_in[row][0] = matchrows[small][id_mat]
            row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]
                nF = False

            elif len(udx)==len(master):
                nF = False


    return master,matchids_in

#
