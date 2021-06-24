import numpy as np
# import matplotlib.pyplot as plt

def matchPSFflc(targname,dir='./',matchtol=3):

    puN = np.genfromtxt(dir+targname+'_allMatchedZPTed_pu.dat',names=True)
    pu = np.genfromtxt(dir+targname+'_allMatchedZPTed_pu.dat')

    psfN = np.genfromtxt(dir + targname + "_PSF2flcTrans_1408.dat",names=True)
    psf = np.genfromtxt(dir + targname + "_PSF2flcTrans_1408.dat")

    id_pu = np.zeros((len(puN),1))
    id_pu[:,0] = np.arange(0,len(puN),1)

    pu = np.hstack((pu,id_pu))
    col_pu = np.array(puN.dtype.names)

    id_psf = np.zeros((len(psfN),1))
    id_psf[:,0] = np.arange(0,len(psfN),1)

    psf = np.hstack((psf,id_psf))
    col_ap = np.array(psfN.dtype.names)

    # Getting columns of x,y of transformed APER to PU
    xA = np.int(np.where(col_ap=='x_Trans')[0])
    yA = np.int(np.where(col_ap=='y_Trans')[0])
    idA = len(col_ap)

    xP = np.int(np.where(col_pu=='xt1_f606w')[0])
    yP = np.int(np.where(col_pu=='yt1_f606w')[0])
    idP = len(col_pu)

    master_in = pu[:,[xP,yP,idP]]
    x,y,idP_mas = 0,1,2

    cat = psf
    matchids = np.zeros((len(master_in),1))

    # nF_out = True
    # matchtol=matchtol

    len_pu = len(pu)
    len_ap = len(psf)

    minLen = np.min([len_pu,len_ap])

    # while nF_out:

    master, matchids = matchlistID(master_in,cat,matchtol,x,y,xA,yA,idA)

        # if len(master)>=int(0.75*minLen):
        #     nF_out = False
        #     print('Minimum Number Reached: %d' % len(master),targname)
        #
        # else:
        #     print('Need More Stars')
        #     print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
        #     matchtol += 1
        #     if matchtol <= 10:
        #         master_in = pu[:,[xP,yP,idP]]
        #         matchids = np.zeros((len(master_in),1))
        #     else:
        #         print("Sacrificing number of stars for quality of matches.")
        #         nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen)

    xP_mas, yP_mas, idP_mas, idA_mas = 0, 1, 2, 3


    idCol_pu = master[:,idP_mas]
    idx_pu = np.asarray(idCol_pu,int)
    reg_pu = pu[idx_pu]

    idCol_ap = master[:,idA_mas]
    idx_ap = np.asarray(idCol_ap,int)
    reg_ap = psf[idx_ap]

    outArr = np.hstack((reg_pu,reg_ap))

    s0 = ' '
    header1 = s0.join(col_pu)
    header1 += ' id_PSFmatch '
    header2 = s0.join(col_ap)
    header2 += ' id_FLCmatch'
    #
    header = header1 + header2

    np.savetxt(dir+targname+'_PSF2flcMatch_1408.dat',outArr,header=header)


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
# matchPSFPu('HOROLOGIUM-I',dir='photUtils0820/',matchtol=3)
