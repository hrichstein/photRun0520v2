import numpy as np
import matplotlib.pyplot as plt


def getRef(targname,filt,dir='./',matchtol=50):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    drcN = np.genfromtxt(dir+'drc_useful_'+targname+'.dat',names=True)

    drc = np.genfromtxt(dir+'drc_useful_'+targname+'.dat')

    idD = np.zeros((len(drc),1))
    idD[:,0] = np.arange(0,len(drc),1)

    drc = np.hstack((drc,idD))

    # flcN = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat',names=True)
    #
    # flc = np.genfromtxt(dir+'magSTDcutAll_'+filt+'.dat')

    flcN = np.genfromtxt(dir+'magSTDcutAll_'+filt+'_mDc.dat',names=True)

    flc = np.genfromtxt(dir+'magSTDcutAll_'+filt+'_mDc.dat')

    idF = np.zeros((len(flc),1))
    idF[:,0] = np.arange(0,len(flc),1)

    flc = np.hstack((flc,idF))

    colDs = np.array(drcN.dtype.names)
    colFs = np.array(flcN.dtype.names)

    idFc = len(colFs) # which column index the id would be
    idDc = len(colDs) # which column index the id would be

    magStr = 'magr'+fils
    xStr = 'x'+fils
    yStr = 'y'+fils

    # Picking  30 brightest stars
    d30 = np.argsort(drcN[magStr])[:30]
    drc30 = drc[d30]

    f30 = np.argsort(flcN['mean'])[:30]
    flc30 = flc[f30]

    # To know which columns have which info
    # Moving from xt,yt space in FLCs to x,y space in DRCs
    xF = np.int(np.where(colFs=='xt1')[0])
    yF = np.int(np.where(colFs=='yt1')[0])
    magF = np.int(np.where(colFs=='mean')[0])
    stdF = np.int(np.where(colFs=='stdev')[0])

    # Getting column of x,y in appropriate filter for DRCs
    xD = np.int(np.where(colDs==xStr)[0])
    yD = np.int(np.where(colDs==yStr)[0])
    magD = np.int(np.where(colDs==magStr)[0])

    master_in = flc30[:,[xF,yF,magF,stdF,idFc]]
    x,y,magrF,stdF_mas,idF_mas = 0,1,2,3,4

    cat = drc30
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    matchtol=matchtol
    while nF_out:
        # master, matchids = matchlistID(master_in,match_arr,matchtol,x1,y1,x2,y2,id_mat)
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrF,stdF_mas,xD,yD,magD,idDc)

        # x,y,magrF,stdF_mas,xD,yD,magD,idDc are all indices, not actual values

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname,filt)
        else:
            print('Need More Stars')
            master_in = flc30[:,[xF,yF,magF,stdF,idFc]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 50

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
    header = 'xF yF magrF stdF idF idD xD yD magrD'
    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f'

    outName = dir+'flcDRCref_'+filt

    # np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    np.savetxt(outName+'_mDc.dat',outArr,header=header,fmt=form)

    # Plotting Section

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xo],outArr[:,yo],label='FLC',s=50)
    ax.scatter(outArr[:,xD],outArr[:,yD],label='DRC',s= 30)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.savefig(outName+'_mDc.png',dpi=600,bbox_inches='tight')

    plt.close()


    return None


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


#
