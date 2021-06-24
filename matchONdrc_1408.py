import numpy as np
# import matplotlib.pyplot as plt

def matchONdrc(targname,dir='./',matchtol=3):

    # drc_dir  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
    # oldN = np.genfromtxt(drc_dir + 'HOROLOGIUM-I_sfErr.dat',names=True)
    # old = np.genfromtxt(drc_dir + 'HOROLOGIUM-I_sfErr.dat')

    oldN = np.genfromtxt(dir + targname + "_oldDRC2newTrans_1408.dat",names=True)
    old = np.genfromtxt(dir + targname + "_oldDRC2newTrans_1408.dat")

    newN = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat',names=True)
    new = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat')

    id_old = np.zeros((len(oldN),1))
    id_old[:,0] = np.arange(0,len(oldN),1)

    old = np.hstack((old,id_old))
    col_old = np.array(oldN.dtype.names)

    id_new = np.zeros((len(newN),1))
    id_new[:,0] = np.arange(0,len(newN),1)

    new = np.hstack((new,id_new))
    col_new = np.array(newN.dtype.names)

    # Getting columns of x,y of transformed OLD to NEW
    xO = np.int(np.where(col_old=='x_Trans')[0])
    yO = np.int(np.where(col_old=='y_Trans')[0])
    idO = len(col_old)

    xN = np.int(np.where(col_new=='xcenter_f606w')[0])
    yN = np.int(np.where(col_new=='ycenter_f606w')[0])
    idN = len(col_new)

    master_in = old[:,[xO,yO,idO]]
    x,y,idO_mas = 0,1,2

    cat = new
    matchids = np.zeros((len(master_in),1))

    # nF_out = True
    # matchtol=matchtol

    len_old = len(old)
    len_new = len(new)

    minLen = np.min([len_old,len_new])

    # while nF_out:

    master, matchids = matchlistID(master_in,cat,matchtol,x,y,xN,yN,idN)

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

    xO_mas, yO_mas, idO_mas, idN_mas = 0, 1, 2, 3


    idCol_old = master[:,idO_mas]
    idx_old = np.asarray(idCol_old,int)
    reg_old = old[idx_old]

    idCol_new = master[:,idN_mas]
    idx_new = np.asarray(idCol_new,int)
    reg_new = new[idx_new]

    outArr = np.hstack((reg_old,reg_new))

    s0 = ' '
    header1 = s0.join(col_old)
    header1 += ' id_DRCmatchN '
    header2 = s0.join(col_new)
    header2 += ' id_DRCmatchO'
    #
    header = header1 + header2

    np.savetxt(dir+targname+'_old2newDRCMatch_1408.dat',outArr,header=header)


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
# matchAPERPu('HOROLOGIUM-I',dir='photUtils0820/',matchtol=3)
