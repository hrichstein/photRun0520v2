import numpy as np

file = np.genfromtxt('../sfdCorrection.dat',names=True)
colNs = np.array(file.dtype.names)

fileS = np.genfromtxt('../sfdCorrection.dat',dtype=str)

targ = np.int(np.where(colNs=='TARGNAME')[0])
f1 = np.int(np.where(colNs=='FILTER1')[0])
f2 = np.int(np.where(colNs=='FILTER2')[0])

targnames= fileS[:,targ]
filt1 = fileS[:,f1]
filt2 = fileS[:,f2]

distFile = np.loadtxt('../nameslist.txt',delimiter=',',usecols=0,dtype=str)
distFf = np.loadtxt('../nameslist.txt',delimiter=',',usecols=5)
distNames = distFile
distDM = distFf

def applyRedDm(targname,dir='./'):

    zptFile = np.genfromtxt(dir+targname+'_fullCat.dat',names=True)
    zptCat = np.genfromtxt(dir+targname+'_fullCat.dat')
    colZs = np.array(zptFile.dtype.names)

    vMag = zptFile['magr_f606w']
    iMag = zptFile['magr_f814w']

    vCorr = file['V_EBV']
    iCorr = file['I_EBV']

    targIdx = []
    # print(targnames)
    for cc,targ in enumerate(targnames):
        if targ==targname:
            targIdx.append(cc)

    for dd,dname in enumerate(distNames):
        if dname==targname:
            distMod=distDM[dd]

    # print(targIdx)

    # Finding which 2 lines correspond to this target (one for each filter)
    # targIdx1 = np.int(np.where(targnames==targname)[0])
    # targIdx2 = np.int(np.where(targnames==targname)[1])
    # # getting the filter information for both of these lines
    tempArr1 = filt1[targIdx]
    tempArr2 = filt2[targIdx]
    # # getting the two corrections for each line
    tempV = vCorr[targIdx]
    tempI = iCorr[targIdx]
    #
    if (tempArr1[0]=='F606W') or (tempArr2[0]=='F606W'):
        v_corr = tempV[0]
        i_corr = tempI[1]
    elif (tempArr1[0]=='F814W') or (tempArr2[0]=='F814W'):
        i_corr = tempI[0]
        v_corr = tempV[1]

    newCols = np.zeros((len(zptFile),4))

    newCols[:,0] = vMag-v_corr
    newCols[:,1] = iMag-i_corr
    newCols[:,2] = vMag-v_corr-distMod
    newCols[:,3] = iMag-i_corr-distMod

    outArr = np.hstack((zptCat,newCols))

    s0 = ' '
    header = s0.join(colZs)
    header += ' magRed_f606w magRed_f814w magAbs_f606w magAbs_f814w'

    np.savetxt(dir+targname+'_absMag.dat',outArr,header=header)


    return None
