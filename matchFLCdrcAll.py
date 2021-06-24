import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits

from linear6d import *
from linTiter_pu import *
from getMatchedFLCdrc0707 import *

def whichIter(targname,filt='F606W',dir='./',matchtol=10):

    diag = open(dir+'iterDiags_all.dat','w')

    iter = int(0)
    kI = True

    nS_list = []
    mOff_list = []
    it_list = []

    while kI:

        numStars, meanOffset = getRef_i(targname,filt=filt,dir=dir,matchtol=matchtol,stdTol=2.5,iter=iter)

        diag.write('Iter:' + str(iter) + '\n')
        diag.write('Number Ref Stars:' + str(numStars) + '\n')
        diag.write('Median Offset:' + str(meanOffset) + '\n')

        nS_list = np.append(nS_list,[numStars])
        mOff_list = np.append(mOff_list,[meanOffset])
        it_list = np.append(it_list,[iter])

        if iter>=int(3): # may need to change to 4 in some cases?? (ex. Hydra2)
            maxStars = it_list[np.argsort(nS_list)[-1]] # iteration that had the largest number of reference stars found, meaning the PREVIOUS iteration's linTrans was best
            minDiff = it_list[np.argsort(mOff_list)[0]]

            print('MaxStars:',maxStars)
            print('minDiff:',minDiff)

            if maxStars==minDiff: #index position
                takeRun = it_list[int(minDiff)-1] # iteration to take
                kI = False
            else:
                if mOff_list[np.argsort(mOff_list)[0]] <= 0.5:
                    takeRun = it_list[int(minDiff)-1]
                    kI = False
            # elif:
            #     np.logical_and(2 >= (it_list[np.argsort(nS_list)[-1]] - it_list[np.argsort(nS_list)[-2]]),minDiff)


        if kI!= False:
            linFLC2drc(targname,filt=filt,dir=dir,iter=iter)
        iter += 1
        if matchtol >= 5:
            matchtol = matchtol - 2
        if iter>=10:
            kI=False
            takeRun = it_list[int(minDiff)-1]
            print('Need to check the diagnostic file.')

    file_str = dir+filt+"_ALL_drcTrans_"+str(int(takeRun))+".dat"

    diag.write('Took Run: ' + str(takeRun) + '\n')
    diag.write('Output in: ' + file_str)
    diag.close()
    # return file_str, nS_list, mOff_list, it_list, out
    return file_str


def doIterMatchDRC(targname,filt='F606W',dir='./',matchtol=2.5,stdTol=5):

    match_file = whichIter(targname,filt=filt,dir=dir)

    print(match_file)

    getMatch(targname,match_file,filt=filt,dir=dir,matchtol=matchtol,stdTol=stdTol)

    return None

def getMatch(targname,file,filt='F606W',dir='./',matchtol=2.5,stdTol=5):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    magdStr = 'magr'+fils
    xdStr = 'xcenter'+fils
    ydStr = 'ycenter'+fils

    magfStr = 'magZPT'+fils
    stdev = 'stdev'+fils

    drcN = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat',names=True)
    drc = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat')

    idD = np.zeros((len(drc),1))
    idD[:,0] = np.arange(0,len(drc),1)

    drc = np.hstack((drc,idD))
    colDs = np.array(drcN.dtype.names)

    # Getting column of x,y in appropriate filter for DRCs
    xD = np.int(np.where(colDs==xdStr)[0])
    yD = np.int(np.where(colDs==ydStr)[0])
    magD = np.int(np.where(colDs==magdStr)[0])

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
    magF = np.int(np.where(colFs==magfStr)[0])
    stdF = np.int(np.where(colFs==stdev)[0])

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

        if len(master)>=int(0.2*len(cat)):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname,filt)
        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            master_in = flc[:,[xF,yF,magF,stdF,idFc]]
            matchids = np.zeros((len(master_in),1))
            matchtol += 2.5

        if matchtol >= 20:
            print("Sacrificing number of stars for quality of matches.")
            nF_out = False

    master = np.hstack((master,matchids))
    print(targname, filt, len(master))

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

    outName = dir+'flcDRCmatch_ALL_'+filt

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

    # Outputting the relevant information columns

    idCol_flc = outArr[:,idF]
    idxFLC = np.asarray(idCol_flc,int)
    regFLC = flc[idxFLC]

    idCol_drc = outArr[:,idD]
    idxDRC = np.asarray(idCol_drc,int)
    regDRC = drc[idxDRC]

    out2 = np.hstack((regFLC,regDRC))

    s0 = ' '
    headerF = s0.join(colFs)
    headerF += ' idf_FDmatch '

    headerD = s0.join(colDs)
    headerD += ' idd_FDmatch'

    header = headerF+headerD

    np.savetxt(dir+targname+'_matchedFLC_DRC_on_'+filt+'.dat',out2,header=header)

    return None
#


def getRef_i(targname,filt='F606W',dir='./',matchtol=5,stdTol=2.5,iter=0):

    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'

    drcN = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat',names=True)
    drc = np.genfromtxt(dir+targname+'_filtMatchDRC_pU.dat')

    idD = np.zeros((len(drc),1))
    idD[:,0] = np.arange(0,len(drc),1)

    drc = np.hstack((drc,idD))
    colDs = np.array(drcN.dtype.names)

    if iter==0:
        flcN = np.genfromtxt(dir+targname+"_allMatchedZPTed_pu.dat",names=True)
        flc = np.genfromtxt(dir+targname+"_allMatchedZPTed_pu.dat")

        colFs = np.array(flcN.dtype.names)

        xfStr = 'xt1'+fils
        yfStr = 'yt1'+fils

        xF = np.int(np.where(colFs==xfStr)[0])
        yF = np.int(np.where(colFs==yfStr)[0])

    else:

        flcN = np.genfromtxt(dir+filt+"_ALL_drcTrans_{}.dat".format(iter-1),names=True)
        flc = np.genfromtxt(dir+filt+"_ALL_drcTrans_{}.dat".format(iter-1))

        colFs = np.array(flcN.dtype.names)

        xstr = 'xDRC_'+ str(int(iter-1))
        ystr = 'yDRC_'+ str(int(iter-1))

        xF = np.int(np.where(colFs==xstr)[0])
        yF = np.int(np.where(colFs==ystr)[0])


    idF = np.zeros((len(flc),1))
    idF[:,0] = np.arange(0,len(flc),1)

    magfStr = 'magZPT'+fils
    stdev = 'stdev'+fils

    flc = np.hstack((flc,idF))

    idFc = len(colFs) # which column index the id would be
    idDc = len(colDs) # which column index the id would be

    magdStr = 'magr'+fils
    xdStr = 'xcenter'+fils
    ydStr = 'ycenter'+fils

    # Picking 50 brightest stars
    d50 = np.argsort(drcN[magdStr])[:50]
    drc50 = drc[d50]

    f50 = np.argsort(flcN[magfStr])[:50]
    flc50 = flc[f50]

    # To know which columns have which info

    magF = np.int(np.where(colFs==magfStr)[0])
    stdF = np.int(np.where(colFs==stdev)[0])

    # Getting column of x,y in appropriate filter for DRCs
    xD = np.int(np.where(colDs==xdStr)[0])
    yD = np.int(np.where(colDs==ydStr)[0])
    magD = np.int(np.where(colDs==magdStr)[0])

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
    mean_off = np.median(tot_off)

    header = 'xF yF magrF stdF idF idD xD yD magrD'
    form = '%1.5f %1.5f %1.4f %1.5f %d %d %1.5f %1.5f %1.4f'

    # changed to reflect iteration
    outName = dir+'flcDRCref_ALL_'+str(iter)

    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    # Plotting Section

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xD],outArr[:,yD],label='DRC',s= 20)
    ax.scatter(outArr[:,xo],outArr[:,yo],label='FLC',s=5)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.close()

    # returning number of reference stars and average (x,y) offset
    print(len(outArr), np.round(mean_off,4))

    return len(outArr), mean_off

def linFLC2drc(targname,filt='F606W',dir='./',iter=0):

    # iter = int(iter)
    if filt=='F606W':
        fils = '_f606w'
    elif filt=='F814W':
        fils = '_f814w'


    file = np.genfromtxt(dir+'flcDRCref_ALL_'+str(iter)+'.dat',names=True)
    fileCat = np.genfromtxt(dir+'flcDRCref_ALL_'+str(iter)+'.dat')

    colNs = np.array(file.dtype.names)

    match_arr = np.zeros((len(file),2))
    xt = np.int(np.where(colNs=='xF')[0])
    yt = np.int(np.where(colNs=='yF')[0])

    match_arr[:,0] = fileCat[:,xt]
    match_arr[:,1] = fileCat[:,yt]

    master_arr = np.zeros((len(file),2))
    xt1 = np.int(np.where(colNs=='xD')[0])
    yt1 = np.int(np.where(colNs=='yD')[0])

    master_arr[:,0] = fileCat[:,xt1]
    master_arr[:,1] = fileCat[:,yt1]

    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    if iter==0:

        all = np.genfromtxt(dir+targname+"_allMatchedZPTed_pu.dat",names=True)
        allCat = np.genfromtxt(dir+targname+"_allMatchedZPTed_pu.dat")

        colAs = np.array(all.dtype.names)

        x_bt = np.int(np.where(colAs=='xt1' +fils)[0])
        y_bt = np.int(np.where(colAs=='yt1' +fils)[0])


    else:

        all = np.genfromtxt(dir+filt+"_ALL_drcTrans_"+str(iter-1)+'.dat',names=True)
        allCat = np.genfromtxt(dir+filt+"_ALL_drcTrans_"+str(iter-1)+'.dat')


        colAs = np.array(all.dtype.names)
        itx = 'xDRC_'+ str(iter-1)
        ity = 'yDRC_'+ str(iter-1)

        x_bt = np.int(np.where(colAs==itx)[0])
        y_bt = np.int(np.where(colAs==ity)[0])


    s0=' '
    header = s0.join(colAs)

    all_arr = np.zeros((len(all),2))

    all_arr[:,0] = allCat[:,x_bt]
    all_arr[:,1] = allCat[:,y_bt]

    iterStr = str(iter)
    outName = dir+filt+"_ALL_drcTrans_"+iterStr

    try:
        new_match, new_all = test_linear(match_arr[:,0],match_arr[:,1], master_arr[:,0], master_arr[:,1], weights, weights, all_arr[:,0],all_arr[:,1])

        makePlot(targname,filt,match_arr[:,0],match_arr[:,1],\
        new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],label_1='Original in FLC',label_2='New in FLC2DRC',label_3='Original in DRC',outname=outName+'_matchCheck')

        outArr = np.hstack((allCat,new_all))
        header += ' xDRC_'+str(iter) + ' yDRC_'+str(iter)

        np.savetxt(outName+'.dat',outArr,header=header)

    except RuntimeWarning:
        print('Not good enough.',targname,filt,dd)


    return None


def makePlot(targname,filt,x1,y1,x2,y2,x3,y3,label_1,\
    label_2,label_3,outname):

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(x3,y3,label=label_3,s=70)
    ax.scatter(x1,y1,label=label_1,s=50)
    ax.scatter(x2,y2,label=label_2,s=20)

    ax.legend()
    ax.set_title(targname+'_'+filt)

    plt.savefig(outname+'.png',dpi=600,bbox_inches='tight')
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
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]
                nF = False

            elif len(udx)==len(master):
                nF = False

    return master,matchids_in
