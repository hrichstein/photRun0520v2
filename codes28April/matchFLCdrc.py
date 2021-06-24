import numpy as np

from matchlistID import matchlistID


def matchFLCdrc(targname,flcFile,drcFile,dir='./',matchtol=1):

    flcN = np.genfromtxt(dir+flcFile,names=True)
    flc = np.genfromtxt(dir+flcFile)

    drcN = np.genfromtxt(dir+drcFile,names=True)
    drc = np.genfromtxt(dir+drcFile)

    colFs = np.array(flcN.dtype.names)
    colDs = np.array(drcN.dtype.names)

    xF = np.int(np.where(colFs=='x_DRCtrans')[0])
    yF = np.int(np.where(colFs=='y_DRCtrans')[0])
    f606w_mag = np.int(np.where(colFs=='mean_f606w')[0])
    f814w_mag = np.int(np.where(colFs=='mean_f814w')[0])
    f606w_err = np.int(np.where(colFs=='err_f606w')[0])
    f814w_err = np.int(np.where(colFs=='err_f814w')[0])

    xD = np.int(np.where(colDs=='xcenter_f606w')[0])
    yD = np.int(np.where(colDs=='ycenter_f606w')[0])

    idColF = len(colFs)
    newCol = np.zeros((len(flc),1),dtype=int)
    newCol[:,0] = np.arange(0,len(flc),1)

    flc_id = np.hstack((flc,newCol))

    idColD = len(colDs)
    newCol = np.zeros((len(drc),1),dtype=int)
    newCol[:,0] = np.arange(0,len(drc),1)
    drc_id = np.hstack((drc,newCol))

    master_in = drc_id[:,[idColD,xD,yD]]

    idD, xd, yd = 0, 1, 2

    cat = flc_id

    nF_out = True

    minLen = len(drc_id)

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,xd,yd,xF,yF,
                                       idColF)

        if len(master)>=int(0.65*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,
                                                             len(master)))

            matchtol += 1
            if matchtol <= 4:
                master_in = drc_id[:,[idColD,xD,yD]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen)

    idD, xd, yd, idF = 0, 1, 2, 3

    idColF = master[:,idF]
    idxF = np.asarray(idColF,int)
    regF = flc[idxF]

    newCols = np.zeros((len(regF),4))
    newCols = regF[:,[f606w_mag,f814w_mag,f606w_err,f814w_err]]

    idColD = master[:,idD]
    idxD = np.asarray(idColD,int)
    regD = drc[idxD]

    allOut = np.hstack((regD,regF))
    headerD = ' '.join(colDs)
    headerF = ' '.join(colFs)

    headerAll = headerD + ' ' + headerF
    np.savetxt(dir+targname+'_fullCat.dat',allOut,header=headerAll)

    shortOut = np.hstack((regD,newCols))
    headerShort = headerD + ' meanFLC_f606w meanFLC_f814w err_f606w err_f814w'

    np.savetxt(dir+targname+'_wErr.dat',shortOut,header=headerShort)

    return None

#
