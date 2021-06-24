import numpy as np
import matplotlib.pyplot as plt

def getONRef(targname,filt='F606W',dir='./',matchtol=10):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    oldN = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I/HOROLOGIUM-I_allMatchedZPTed.dat',names=True)
    old = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I/HOROLOGIUM-I_allMatchedZPTed.dat')

    idO = np.zeros((len(old),1))
    idO[:,0] = np.arange(0,len(old),1)

    newN = np.genfromtxt(dir+'HOROLOGIUM-I_allMatchedZPTed_pu.dat',names=True)
    new = np.genfromtxt(dir+'HOROLOGIUM-I_allMatchedZPTed_pu.dat')

    idN = np.zeros((len(new),1))
    idN[:,0] = np.arange(0,len(new),1)

    old = np.hstack((old,idO))
    new = np.hstack((new,idN))

    colOs = np.array(oldN.dtype.names)
    colNs = np.array(newN.dtype.names)

    idOc = len(colOs)
    idNc = len(colNs)

    magStr = 'magZPT'+ fils
    xStr = 'xt1' + fils
    yStr = 'yt1' + fils
    stdev = 'stdev' + fils

    o30 = np.argsort(oldN[magStr])[:30]
    old30 = old[o30]

    n30 = np.argsort(newN[magStr])[:30]
    new30 = new[n30]

    xO = np.int(np.where(colOs==xStr)[0])
    yO = np.int(np.where(colOs==yStr)[0])
    magO = np.int(np.where(colOs==magStr)[0])
    stdO = np.int(np.where(colOs==stdev)[0])

    xN = np.int(np.where(colNs==xStr)[0])
    yN = np.int(np.where(colNs==yStr)[0])
    magN = np.int(np.where(colNs==magStr)[0])
    stdN = np.int(np.where(colNs==stdev)[0])

    master_in = old30[:,[xO,yO,magO,stdO,idOc]]
    x,y,magrO,stdO_mas,idO_mas = 0,1,2,3,4

    cat = new30
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    matchtol=matchtol

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrO,stdO_mas,xN,yN,magN,idNc)

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname,filt)
        else:
            print('Need More Stars')
            master_in = old30[:,[xO,yO,magO,stdO,idOc]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 5

    master = np.hstack((master,matchids))

    xO_mas, yO_mas, magO_mas, stdO_mas, idO_mas, idN_mas = 0, 1, 2, 3, 4, 5

    newCols = np.zeros((len(master),4))

    idxCol = master[:,idN_mas]
    idxN = np.asarray(idxCol,int)
    regN = new[idxN]

    newCols[:,0] = regN[:,xN]
    newCols[:,1] = regN[:,yN]
    newCols[:,2] = regN[:,magN]
    newCols[:,3] = regN[:,stdN]

    outArr = np.hstack((master,newCols))

    xO, yO, magrO, stdO, idO, idN, xN, yN, magrN, stdN = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    header = 'xO yO magrO stdO idO idN xN yN magrN stdN'
    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f %1.5f'

    outName = dir+'oldNewref_'+filt

    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xO],outArr[:,yO],label='Old FLC',s=50)
    ax.scatter(outArr[:,xN],outArr[:,yN],label='New FLC',s= 30)

    ax.legend()
    ax.set_title(targname+'_matchedin_'+filt)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None

def getONmatch(targname,filt='F606W',dir='./',matchtol=3,stdTol=5):

    if filt=='F606W':
        fils = '_f606w'
        fils2 = '_f814w'
    elif filt=='F814W':
        fils = '_f814w'
        fils2 = '_f606w'

    xStr = 'xt1'+fils
    yStr = 'yt1'+fils

    magStr = 'magZPT'+fils
    stdev = 'stdev'+fils

    magStr2 = 'magZPT'+fils2
    stdev2 = 'stdev'+fils2

    oldN = np.genfromtxt(dir+'HOROLOGIUM-I_ONtrans_F606W.dat',names=True)
    oldCat = np.genfromtxt(dir+'HOROLOGIUM-I_ONtrans_F606W.dat')
    colOs = np.array(oldN.dtype.names)

    idO = np.zeros((len(oldN),1))
    idO[:,0] = np.arange(0,len(oldN),1)

    old = np.hstack((oldCat,idO))
    colOs = np.array(oldN.dtype.names)

    # Getting column of x,y in appropriate filter for DRCs
    xO = np.int(np.where(colOs=='xt1_ON')[0])
    yO = np.int(np.where(colOs=='yt1_ON')[0])
    magO = np.int(np.where(colOs==magStr)[0])
    stdO = np.int(np.where(colOs==stdev)[0])

    magO2 = np.int(np.where(colOs==magStr2)[0])
    stdO2 = np.int(np.where(colOs==stdev2)[0])

    new = np.genfromtxt(dir+'HOROLOGIUM-I_allMatchedZPTed_pu.dat',names=True)
    newCat = np.genfromtxt(dir+'HOROLOGIUM-I_allMatchedZPTed_pu.dat')
    colNs = np.array(new.dtype.names)

    xN = np.int(np.where(colNs==xStr)[0])
    yN = np.int(np.where(colNs==yStr)[0])
    magN = np.int(np.where(colNs==magStr)[0])
    stdN = np.int(np.where(colNs==stdev)[0])

    magN2 = np.int(np.where(colNs==magStr2)[0])
    stdN2 = np.int(np.where(colNs==stdev2)[0])

    idN = np.zeros((len(new),1))
    idN[:,0] = np.arange(0,len(new),1)

    new = np.hstack((newCat,idN))

    idOc = len(colOs) # which column index the id would be
    idNc = len(colNs) # which column index the id would be

    master_in = new[:,[xN,yN,magN,stdN,idNc]]
    x,y,magrN,stdN_mas,idN_mas = 0,1,2,3,4

    cat = old
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    matchtol=matchtol

    # while nF_out:
    master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrN,stdN_mas,xO,yO,magO,idOc,stdTol=stdTol)


        # if len(master)>=int(0.2*len(cat)):
        #     nF_out = False
        #     print('Minimum Number Reached: %d' % len(master),targname,filt)
        # else:
        #     print('Need More Stars')
        #     print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
        #     master_in = new[:,[xN,yN,magN,stdN,idNc]]
        #     matchids = np.zeros((len(master_in),1))
        #     matchtol += 2.5
        #
        # if matchtol >= 20:
        #     print("Sacrificing number of stars for quality of matches.")
        #     nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master))

    xN_mas, yN_mas, magN_mas, stdN_mas, idN_mas, idO_mas = 0, 1, 2, 3, 4, 5

    newCols = np.zeros((len(master),8))

    idxCol = master[:,idO_mas]
    idxO = np.asarray(idxCol,int)
    regO = old[idxO]

    newCols[:,0] = regO[:,xO]
    newCols[:,1] = regO[:,yO]
    newCols[:,2] = regO[:,magO]
    newCols[:,3] = regO[:,stdO]
    newCols[:,4] = regO[:,magO2]
    newCols[:,5] = regO[:,stdO2]

    idxCol = master[:,idN_mas]
    idxN = np.asarray(idxCol,int)
    regN = new[idxN]

    newCols[:,6] = regN[:,magN2]
    newCols[:,7] = regN[:,stdN2]

    outArr = np.hstack((master,newCols))

    xN, yN, magrN, stdN, idN, idO, xO, yO, magrO,stdO = 0, 1, 2, 3, 4, 5, 6, 7, 8,9

    if filt=='F606W':
        header = 'xN yN magr_f606wN std_f606wN idN idO xO yO magr_f606wO std_f606wO magr_f814wO std_f814wO magr_f814wN std_f814wN'

    elif filt=='F814W':
        header = 'xN yN magr_f814wN std_f814wN idN idO xO yO magr_f814wO std_f814wO magr_f606wO std_f606wO magr_f606wN std_f606wN'


    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f %1.5f %1.4f %1.5f %1.5f %1.4f'

    outName = dir+'newOLDmatch_FULL_'+filt

    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)


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
