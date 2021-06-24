import numpy as np
import matplotlib.pyplot as plt

def getRefPSFflc(targname,matchtol=3,dir='./'):

    flcN = np.genfromtxt(dir+targname+'_allMatchedZPTed_pu.dat',names=True)
    flc = np.genfromtxt(dir+targname+'_allMatchedZPTed_pu.dat')

    psfN = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I_pastRuns/matchedFLCpsf2506_all.dat',names=True)
    psf = np.genfromtxt('catRawMags1305/catDir_HOROLOGIUM-I_pastRuns/matchedFLCpsf2506_all.dat')

    idF = np.zeros((len(flc),1))
    idF[:,0] = np.arange(0,len(flc),1)

    flc = np.hstack((flc,idF))
    colFs = np.array(flcN.dtype.names)

    idP = np.zeros((len(psf),1))
    idP[:,0] = np.arange(0,len(psf),1)

    psf = np.hstack((psf,idP))
    colPs = np.array(psfN.dtype.names)

    xF = np.int(np.where(colFs=='xt1_f606w')[0])
    yF = np.int(np.where(colFs=='yt1_f606w')[0])
    idF = len(colFs)

    xP = np.int(np.where(colPs=='xt1_f606w')[0])
    yP = np.int(np.where(colPs=='yt1_f606w')[0])
    idP = len(colPs)

    xP_tru = np.int(np.where(colPs=='xPSF')[0])
    yP_tru = np.int(np.where(colPs=='yPSF')[0])

    f30 = np.argsort(flcN['magZPT_f606w'])[:30]
    ff_30 = flc[f30]

    p30 = np.argsort(psfN['magZPT_f606w'])[:30]
    pp_30 = psf[p30]

    master_in = ff_30[:,[xF,yF,idF]]
    x, y, idF_mas = 0, 1, 2

    cat = pp_30
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xP,yP,idP)

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = ff_30[:,[xF,yF,idF]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 1

    master = np.hstack((master,matchids))

    xF, yF, idF, idP = 0, 1, 2, 3

    newCols = np.zeros((len(master),2))

    idxCol = master[:,idP]
    idxP = np.asarray(idxCol,int)
    regP = psf[idxP]

    newCols[:,0] = regP[:,xP_tru]
    newCols[:,1] = regP[:,yP_tru]

    outArr = np.hstack((master,newCols))

    xF, yF = 0, 1

    header = 'x_f606w y_f606w idFLC idPSF xPSF yPSF'

    outName = dir+'flcPSFref_'+targname
    np.savetxt(outName+'_1408.dat',outArr,header=header)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xF],outArr[:,yF],label='FLC',s= 50)
    ax.scatter(regP[:,xP],regP[:,yP],label='PSF',s=20)

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


# getRefPSFflc('HOROLOGIUM-I',matchtol=1,dir='photUtils0820/')
