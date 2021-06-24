import numpy as np
# import matplotlib.pyplot as plt

def matchFilt(targname,dir='./',matchtol=3):

    f606wN = np.genfromtxt(dir+'magZPTedAll_F606W_mDc.dat',names=True)
    f606w = np.genfromtxt(dir+'magZPTedAll_F606W_mDc.dat')

    # f814wN = np.genfromtxt(dir+'magZPTedAll_F814W.dat',names=True)
    # f814w = np.genfromtxt(dir+'magZPTedAll_F814W.dat')
    f814wN = np.genfromtxt(dir+targname + "_filtTrans_mDc.dat",names=True)
    f814w = np.genfromtxt(dir+targname + "_filtTrans_mDc.dat")

    id606 = np.zeros((len(f606w),1))
    id606[:,0] = np.arange(0,len(f606w),1)

    f606w = np.hstack((f606w,id606))
    col606 = np.array(f606wN.dtype.names)

    id814 = np.zeros((len(f814w),1))
    id814[:,0] = np.arange(0,len(f814w),1)

    f814w = np.hstack((f814w,id814))
    col814 = np.array(f814wN.dtype.names)

    # Getting columns of x,y of transformed F814W to F606W
    xI = np.int(np.where(col814=='x_f606wTrans')[0])
    yI = np.int(np.where(col814=='y_f606wTrans')[0])
    magI = np.int(np.where(col814=='magZPT')[0])
    idI = len(col814)

    xV = np.int(np.where(col606=='xt1')[0])
    yV = np.int(np.where(col606=='yt1')[0])
    magV = np.int(np.where(col606=='magZPT')[0])
    idV = len(col606)

    master_in = f606w[:,[xV,yV,magV,idV]]
    x,y,magrV,idV_mas = 0,1,2,3

    cat = f814w
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    len606 = len(f606w)
    len814 = len(f814w)

    minLen = np.min([len606,len814])

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrV,xI,yI,magI,idI)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 1
            if matchtol <= 10:
                master_in = f606w[:,[xV,yV,magV,idV]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen)

    xV_mas, yV_mas, magV_mas, idV_mas, idI_mas = 0, 1, 2, 3, 4

    newCols = np.zeros((len(master),3))

    idCol606 = master[:,idV_mas]
    idx606 = np.asarray(idCol606,int)
    reg606 = f606w[idx606]

    idCol814 = master[:,idI_mas]
    idx814 = np.asarray(idCol814,int)
    reg814 = f814w[idx814]

    outArr = np.hstack((reg606,reg814))
    header = 'RA_f606w DEC_f606w flux_f606w flags_f606w c_star_f606w mag1_f606w mag2_f606w mag3_f606w mag4_f606w ra1_f606w dec1_f606w ra2_f606w dec2_f606w ra3_f606w dec3_f606w ra4_f606w dec4_f606w xr1_f606w yr1_f606w xr2_f606w yr2_f606w xr3_f606w yr3_f606w xr4_f606w yr4_f606w xc1_f606w yc1_f606w xc2_f606w yc2_f606w xc3_f606w yc3_f606w xc4_f606w yc4_f606w xt1_f606w yt1_f606w xt2_f606w yt2_f606w xt3_f606w yt3_f606w xt4_f606w yt4_f606w mean_f606w stdev_f606w cut_flag_f606w idx_cut_f606w num_abv_std_f606w magZPT_f606w magZPTerr_f606w id_f606w RA_f814w DEC_f814w flux_f814w flags_f814w c_star_f814w mag1_f814w mag2_f814w mag3_f814w mag4_f814w ra1_f814w dec1_f814w ra2_f814w dec2_f814w ra3_f814w dec3_f814w ra4_f814w dec4_f814w xr1_f814w yr1_f814w xr2_f814w yr2_f814w xr3_f814w yr3_f814w xr4_f814w yr4_f814w xc1_f814w yc1_f814w xc2_f814w yc2_f814w xc3_f814w yc3_f814w xc4_f814w yc4_f814w xt1_f814w yt1_f814w xt2_f814w yt2_f814w xt3_f814w yt3_f814w xt4_f814w yt4_f814w mean_f814w stdev_f814w cut_flag_f814w idx_cut_f814w num_abv_std_f814w magZPT_f814w magZPTerr_f814w id_f814w x_f606wTrans y_f606wTrans'

    # np.savetxt(dir+targname+'_allMatchedZPTed.dat',outArr,header=header)
    np.savetxt(dir+targname+'_allMatchedZPTed_mDc.dat',outArr,header=header)


    return None

def matchlistID(master,cat,matchtol,x1,y1,mag1,\
    x2,y2,mag2,id_mat):

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
