import numpy as np
import os
from astropy.io import fits

from refIter_pu import *
from linTiter_pu import *
from getMatchedFLCdrc0707 import *

def whichIter(targname,filt,dir='./',matchtol=10):

    diag = open(dir+'iterDiags_'+filt+'.dat','w')

    iter = int(1)
    kI = True

    nS_list = []
    mOff_list = []
    it_list = []

    takeRun = int(0)

    while kI:

        numStars, meanOffset = getRef_i(targname,filt,dir=dir,matchtol=matchtol,stdTol=2.5,iter=iter)

        diag.write('Iter:' + str(iter) + '\n')
        diag.write('Number Ref Stars:' + str(numStars) + '\n')
        diag.write('Median Offset:' + str(meanOffset) + '\n')

        nS_list = np.append(nS_list,[numStars])
        mOff_list = np.append(mOff_list,[meanOffset])
        it_list = np.append(it_list,[iter])

        if iter>=int(3): # may need to change to 4 in some cases?? (ex. Hydra2)
            maxStars = it_list[np.argsort(nS_list)[-1]] # iteration that had the largest number of reference stars found, meaning the PREVIOUS iteration's linTrans was best
            minDiff = it_list[np.argsort(mOff_list)[0]]

            if maxStars==minDiff:
                takeRun = int(minDiff-1) # iteration to take
                kI = False
            elif np.argsort(mOff_list)[0] <= 1:
                takeRun = int(minDiff-1)
                kI = False

            # if np.logical_and((nS_list[-1] >= nS_list[-2]),(mOff_list[-1] <= mOff_list[-2])) :
            #     kI = False
            #     takeRun = iter
            #     # out = 'Case 1'
            #
            # elif np.logical_and((nS_list[-1] >= nS_list[-2]),(mOff_list[-1] >= mOff_list[-2])) :
            #     kI = False
            #     takeRun = iter-1
            #     # out = 'Case 2'
            #
            # elif np.logical_and((nS_list[-1] <= nS_list[-2]),(mOff_list[-1] >= mOff_list[-2])) :
            #     kI = False
            #     takeRun = iter-1
            #     # out = 'Case 3'

        if kI!= False:
            linFLC2drc_i(targname,filt,dir=dir,iter=iter)
        iter += 1
        if matchtol > 5:
            matchtol = matchtol - 2

        if iter>=10:
            takeRun = int(minDiff-1)
            kI = False

    if takeRun==0:
        file_str = dir+targname+"_"+filt+"_drcTrans_pu.dat"
    # elif takeRun == 1:
    #     file_str = dir+targname+"_"+filt+"_drcTrans.dat"
    else:
        # iterStr = str(iter) # I think I don't have to subtract one here... I hope
        iterStr = str(takeRun)
        file_str = dir+targname+"_"+filt+"_drcTrans"+iterStr+".dat"

    diag.write('Took Run:' + str(takeRun) + '\n')
    diag.write('Output in:' + file_str)
    diag.close()
    # return file_str, nS_list, mOff_list, it_list, out
    return file_str


def doIterMatch(targname,filt,dir='./',matchtol=2.5,stdTol=5):

    match_file = whichIter(targname,filt,dir=dir)

    print(match_file)

    getMatch(targname,filt,match_file,dir=dir,matchtol=matchtol,stdTol=stdTol)

    return None

def getMatch(targname,filt,file,dir='./',matchtol=2.5,stdTol=5):

    if filt=='F606W':
        fils = '_f606w'
        fils2 = '_f814w'
    elif filt=='F814W':
        fils = '_f814w'
        fils2 = '_f606w'

    magStr = 'magr'+fils
    xStr = 'xcenter'+fils
    yStr = 'ycenter'+fils

    magStr2 = 'magr'+fils2

    drcDir = './photUtils20Aug/catDir_'+targname+'/'

    drcN = np.genfromtxt(drcDir+targname+'_filtMatchDRC_pU.dat',names=True)
    drc = np.genfromtxt(drcDir+targname+'_filtMatchDRC_pU.dat')

    idD = np.zeros((len(drc),1))
    idD[:,0] = np.arange(0,len(drc),1)

    drc = np.hstack((drc,idD))
    colDs = np.array(drcN.dtype.names)

    # Getting column of x,y in appropriate filter for DRCs
    xD = np.int(np.where(colDs==xStr)[0])
    yD = np.int(np.where(colDs==yStr)[0])
    magD = np.int(np.where(colDs==magStr)[0])
    magD2 = np.int(np.where(colDs==magStr2)[0])

    # Things affected by iter changing
    # flcN = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter)+".dat",names=True)
    # flc = np.genfromtxt(dir+targname+"_"+filt+"_drcTrans"+str(iter)+".dat")
    # file name is of the form dir+targname+"_"+filt+"_drcTrans"+iterStr+".dat"

    flcN = np.genfromtxt(file,names=True)
    flc = np.genfromtxt(file)

    colFs = np.array(flcN.dtype.names)

    xstr = colFs[-2]
    ystr = colFs[-1]

    xF = np.int(np.where(colFs==xstr)[0])
    yF = np.int(np.where(colFs==ystr)[0])
    magF = np.int(np.where(colFs=='mean')[0])
    stdF = np.int(np.where(colFs=='stdev')[0])

    idF = np.zeros((len(flc),1))
    idF[:,0] = np.arange(0,len(flc),1)

    flc = np.hstack((flc,idF))

    idFc = len(colFs) # which column index the id would be
    idDc = len(colDs) # which column index the id would be

    master_in = flc[:,[xF,yF,magF,stdF,idFc]]
    x,y,magrF,stdF_mas,idF_mas = 0,1,2,3,4

    cat = drc
    matchids = np.zeros((len(master_in),1))

    nF_out = True

    matchtol=matchtol

    while nF_out:
        # master, matchids = matchlistID(master_in,match_arr,matchtol,x1,y1,x2,y2,id_mat)
        master, matchids = matchlistID(master_in,cat,matchtol,x,y,magrF,stdF_mas,xD,yD,magD,idDc,stdTol=stdTol)

        # x,y,magrF,stdF_mas,xD,yD,magD,idDc are all indices, not actual values

        if len(master)>=int(0.1*len(cat)):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname,filt)
        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 2.5
            if matchtol <= 20:
                master_in = flc[:,[xF,yF,magF,stdF,idFc]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF = False

        # if matchtol >= 20:
        #     print("Sacrificing number of stars for quality of matches.")
        #     nF_out = False

    master = np.hstack((master,matchids))
    print(targname, filt, len(master))

    xF_mas, yF_mas, magF_mas, stdF_mas, idF_mas, idD_mas = 0, 1, 2, 3, 4, 5

    newCols = np.zeros((len(master),4))

    idxCol = master[:,idD_mas]
    idxD = np.asarray(idxCol,int)
    regD = drc[idxD]

    newCols[:,0] = regD[:,xD]
    newCols[:,1] = regD[:,yD]
    newCols[:,2] = regD[:,magD]
    newCols[:,3] = regD[:,magD2]

    outArr = np.hstack((master,newCols))

    xo, yo, magrF, stdF, idF, idD, xD, yD, magrD = 0, 1, 2, 3, 4, 5, 6, 7, 8
    header = 'xF yF magrF stdF idF idD xD yD '
    header += 'magrD'+fils
    header += ' magrD'+fils2
    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f %1.4f'

    outName = dir+targname+'_flcDRCmatch_'+filt

    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    # Plotting Section

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xD],outArr[:,yD],label='DRC',s= 50)
    ax.scatter(outArr[:,xo],outArr[:,yo],label='FLC',s=20)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.close()

    fig, ax = plt.subplots(figsize=(6,4))

    ax.scatter(outArr[:,magrD],outArr[:,magrD]-outArr[:,magrF],s=20)
    ax.set_xlabel('DRC mag')
    ax.set_ylabel('DRC - Mean mag')

    ax.set_title(targname+'_'+filt)

    outName = dir+'magDiffs_'+filt

    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

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

            matchids_in[row][0] = matchrows[small][id_mat]
            row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            print('Tripping UDX')
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]
                print(len(master),len(matchids_in))
                nF = False

            elif (len(udx)==len(master)) and (len(udx)==len(matchids_in)) :
                nF = False

    return master,matchids_in

#
