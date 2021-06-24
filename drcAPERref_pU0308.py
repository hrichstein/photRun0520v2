import numpy as np
import matplotlib.pyplot as plt

def getRefApPu(targname,matchtol=3,dir='photUtils0820/'):

    pUN = np.genfromtxt(dir+'HOROLOGIUM-I_filtMatch_pU.dat',names=True)
    pU_cat = np.genfromtxt(dir+'HOROLOGIUM-I_filtMatch_pU.dat')

    # dir2 = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
    # aperN = np.genfromtxt(dir2+'HOROLOGIUM_CF.2.TOSEND.CAT',names=True)
    # aper_cat = np.genfromtxt(dir2+'HOROLOGIUM_CF.2.TOSEND.CAT')
    dir2 = 'catRawMags1305/catDir_HOROLOGIUM-I_pastRuns/'
    aperN = np.genfromtxt(dir2+'matchedFLCpsf2506_all.dat',names=True)
    aper_cat = np.genfromtxt(dir2+'matchedFLCpsf2506_all.dat')

    id_pu = np.zeros((len(pUN),1))
    id_pu[:,0] = np.arange(0,len(pUN),1)

    pU_cat = np.hstack((pU_cat,id_pu))
    puNames = np.array(pUN.dtype.names)

    id_ap = np.zeros((len(aperN),1))
    id_ap[:,0] = np.arange(0,len(aperN),1)

    aper_cat = np.hstack((aper_cat,id_ap))
    apNames = np.array(aperN.dtype.names)

    xP = np.int(np.where(puNames=='xcenter_f606w')[0])
    yP = np.int(np.where(puNames=='ycenter_f606w')[0])
    idP = len(puNames)

    xA = np.int(np.where(apNames=='xDRC_mat_f606w')[0])
    yA = np.int(np.where(apNames=='yDRC_mat_f606w')[0])
    idA = len(apNames)

    xA_tru = np.int(np.where(apNames=='xAPER')[0])
    yA_tru = np.int(np.where(apNames=='yAPER')[0])

    p50 = np.argsort(pUN['magr_f606w'])[:50]
    fp_50 = pU_cat[p50]

    a50 = np.argsort(aperN['m606cAPER'])[:50]
    fa_50 = aper_cat[a50]

    master_in = fp_50[:,[xP,yP,idP]]
    x, y, idP_mas = 0, 1, 2

    cat = fa_50
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xA,yA,idA)

        if len(master)>=int(6): # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = fp_50[:,[xP,yP,idP]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 1

    master = np.hstack((master,matchids))

    xPU, yPU, idPU, idAP = 0, 1, 2, 3

    newCols = np.zeros((len(master),2))

    idxCol = master[:,idAP]
    idxAP = np.asarray(idxCol,int)
    regAP = aper_cat[idxAP]

    newCols[:,0] = regAP[:,xA_tru]
    newCols[:,1] = regAP[:,yA_tru]

    outArr = np.hstack((master,newCols))

    xP, yP, idP, idA, xA_mas, yA_mas = 0, 1, 2, 3, 4, 5

    header = 'xPU_f606w yPU_f606w idPU idAP xAPER yAPER'

    outName = dir+'drcAPERref_'+targname
    np.savetxt(outName+'_pU.dat',outArr,header=header)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xP],outArr[:,yP],label='PhotU',s= 50)
    ax.scatter(regAP[:,xA],regAP[:,yA],label='APER',s=20)

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outName+'_pU.png',dpi=600,bbox_inches='tight')

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


getRefApPu('HOROLOGIUM-I',matchtol=1,dir='photUtils0820/')
