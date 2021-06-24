import numpy as np
import matplotlib.pyplot as plt


def getRef_i(targname,filt,dir='./',matchtol=5,stdTol=2.5,iter=1):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    drcN = np.genfromtxt(dir+'drc_useful_'+targname+'.dat',names=True)
    drc = np.genfromtxt(dir+'drc_useful_'+targname+'.dat')

    idD = np.zeros((len(drc),1))
    idD[:,0] = np.arange(0,len(drc),1)

    drc = np.hstack((drc,idD))
    colDs = np.array(drcN.dtype.names)

    # if iter==1:
    #     flcN = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans.dat",names=True)
    #     flc = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans.dat")
    #
    #     colFs = np.array(flcN.dtype.names)
    #
    #     xF = np.int(np.where(colFs=='xDRC')[0])
    #     yF = np.int(np.where(colFs=='yDRC')[0])
    # else:
    #     flcN = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans{0:d}.dat".format(iter-1),names=True)
    #     flc = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans{0:d}.dat".format(iter-1))
    #
    #     colFs = np.array(flcN.dtype.names)
    #
    #     xstr = 'xDRC_'+ str(int(iter-1))
    #     ystr = 'yDRC_'+ str(int(iter-1))
    #
    #     xF = np.int(np.where(colFs==xstr)[0])
    #     yF = np.int(np.where(colFs==ystr)[0])
    if iter==1:
        flcN = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans_mDc.dat",names=True)
        flc = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans_mDc.dat")

        colFs = np.array(flcN.dtype.names)

        xF = np.int(np.where(colFs=='xDRC')[0])
        yF = np.int(np.where(colFs=='yDRC')[0])
    else:
        flcN = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans{0:d}_mDc.dat".format(iter-1),names=True)
        flc = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans{0:d}_mDc.dat".format(iter-1))

        colFs = np.array(flcN.dtype.names)

        xstr = 'xDRC_'+ str(int(iter-1))
        ystr = 'yDRC_'+ str(int(iter-1))

        xF = np.int(np.where(colFs==xstr)[0])
        yF = np.int(np.where(colFs==ystr)[0])


    # These line differ from the original, which had xt1, etc., since this is the second iteration of this

        ######

    idF = np.zeros((len(flc),1))
    idF[:,0] = np.arange(0,len(flc),1)

    flc = np.hstack((flc,idF))

    idFc = len(colFs) # which column index the id would be
    idDc = len(colDs) # which column index the id would be

    magStr = 'magr'+fils
    xStr = 'x'+fils
    yStr = 'y'+fils

    # Picking 50 brightest stars
    d50 = np.argsort(drcN[magStr])[:50]
    drc50 = drc[d50]

    f50 = np.argsort(flcN['mean'])[:50]
    flc50 = flc[f50]

    # To know which columns have which info

    magF = np.int(np.where(colFs=='mean')[0])
    stdF = np.int(np.where(colFs=='stdev')[0])

    # Getting column of x,y in appropriate filter for DRCs
    xD = np.int(np.where(colDs==xStr)[0])
    yD = np.int(np.where(colDs==yStr)[0])
    magD = np.int(np.where(colDs==magStr)[0])

    master_in = flc50[:,[xF,yF,magF,stdF,idFc]]
    x,y,magrF,stdF_mas,idF_mas = 0,1,2,3,4

    cat = drc50
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    matchtol=matchtol
    while nF_out:
        # master, matchids = matchlistID(master_in,match_arr,matchtol,x1,y1,x2,y2,id_mat)
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrF,stdF_mas,xD,yD,magD,idDc,stdTol=stdTol)

        # x,y,magrF,stdF_mas,xD,yD,magD,idDc are all indices, not actual values

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname,filt)
        else:
            print('Need More Stars')
            master_in = flc50[:,[xF,yF,magF,stdF,idFc]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 5

    master = np.hstack((master,matchids))

    xF_mas, yF_mas, magF_mas, stdF_mas, idF_mas, idD_mas = 0, 1, 2, 3, 4, 5

    newCols = np.zeros((len(master),3))

    idxCol = master[:,idD_mas]
    idxD = np.asarray(idxCol,int)
    regD = drc[idxD]

    newCols[:,0] = regD[:,xD]
    newCols[:,1] = regD[:,yD]
    newCols[:,2] = regD[:,magD]

    outArr = np.hstack((master,newCols))

    xo, yo, magrF, stdF, idF, idD, xD, yD, magrD = 0, 1, 2, 3, 4, 5, 6, 7, 8

    xoff = outArr[:,xo] - outArr[:,xD]
    yoff = outArr[:,yo] - outArr[:,yD]
    tot_off = np.sqrt(xoff**2 + yoff**2)
    mean_off = np.mean(tot_off)

    header = 'xF yF magrF stdF idF idD xD yD magrD'
    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f'

    # changed to reflect iteration
    outName = dir+'flcDRCref_'+filt+'_'+str(iter)

    # np.savetxt(outName+'.dat',outArr,header=header,fmt=form)
    np.savetxt(outName+'_mDc.dat',outArr,header=header,fmt=form)

    # Plotting Section

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xD],outArr[:,yD],label='DRC',s= 20)
    ax.scatter(outArr[:,xo],outArr[:,yo],label='FLC',s=5)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outName+'_mDc.png',dpi=600,bbox_inches='tight')

    plt.close()

    # returning number of reference stars and average (x,y) offset
    print(len(outArr), np.round(mean_off,4))

    return len(outArr), mean_off


def matchlistID(master,cat,matchtol,x1,y1,mag1,std1,\
    x2,y2,mag2,id_mat,stdTol=3):

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

            # Magnitude checking
            if abs(matchrows[small][mag2] - master[row][mag1]) < (stdTol*master[row][std1]):
                matchids_in[row][0] = matchrows[small][id_mat]
                row += 1
            else:
                master = np.delete(master,row,0)
                matchids_in = np.delete(matchids_in,row,0)
            # Failed terribly!!
            # magDif = np.zeros((len(matchrows),1))
            # for mm in range(len(matchrows)):
            #     magDif[mm] = master[row][magrF] - matchrows[mm][magD]
            #     small = np.argmin(magDif)
            #     matchids[row][0] = matchrows[small][idDc]
            # row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):

                # print(len(udx),len(master))
                master = master[udx]
                matchids_in = matchids_in[udx]

                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(matchtol,len(master)))
                nF = False

            elif len(udx)==len(master):
                # print(len(udx),len(master))
                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(matchtol,len(master)))
                nF = False

    return master,matchids_in
