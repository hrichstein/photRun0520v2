import numpy as np
import matplotlib.pyplot as plt

def getRefdrcON(targname,matchtol=3,dir='./'):

    drc_dir  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
    drcON = np.genfromtxt(drc_dir + 'HOROLOGIUM-I_sfErr.dat',names=True)
    drcOld = np.genfromtxt(drc_dir + 'HOROLOGIUM-I_sfErr.dat')

    drcN = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat',names=True)
    drcNew = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat')

    idO = np.zeros((len(drcOld),1))
    idO[:,0] = np.arange(0,len(drcOld),1)

    drcOld = np.hstack((drcOld,idO))
    colOs = np.array(drcON.dtype.names)

    idN = np.zeros((len(drcNew),1))
    idN[:,0] = np.arange(0,len(drcNew),1)

    drcNew = np.hstack((drcNew,idN))
    colNs = np.array(drcN.dtype.names)

    xO = np.int(np.where(colOs=='x_v')[0])
    yO = np.int(np.where(colOs=='y_v')[0])
    idO = len(colOs)

    xN = np.int(np.where(colNs=='xcenter_f606w')[0])
    yN = np.int(np.where(colNs=='ycenter_f606w')[0])
    idN = len(colNs)

    o30 = np.argsort(drcON['magRaw_v'])[:30]
    oo_30 = drcOld[o30]

    n30 = np.argsort(drcN['magr_f606w'])[:30]
    nn_30 = drcNew[n30]

    master_in = nn_30[:,[xN,yN,idN]]
    x, y, idN_mas = 0, 1, 2

    cat = oo_30
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xO,yO,idO)

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = nn_30[:,[xN,yN,idN]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 1

    master = np.hstack((master,matchids))

    xN, yN, idN, idO = 0, 1, 2, 3

    newCols = np.zeros((len(master),2))

    idxCol = master[:,idO]
    idxO = np.asarray(idxCol,int)
    regO = drcOld[idxO]

    newCols[:,0] = regO[:,xO]
    newCols[:,1] = regO[:,yO]

    outArr = np.hstack((master,newCols))

    xN, yN = 0, 1

    header = 'x_new y_new idNew idOld x_old y_old'

    outName = dir+'oldNewDRCref_'+targname
    np.savetxt(outName+'_1408.dat',outArr,header=header)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xN],outArr[:,yN],label='New DRC',s= 50)
    ax.scatter(regO[:,xO],regO[:,yO],label='Old DRC',s=20)

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outName+'_1408.png',dpi=600,bbox_inches='tight')

    plt.close()

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
