import numpy as np

def match4cStar606(targname,dir='./',matchtol=3):

    oldN = np.genfromtxt('catRawMags1305/catDir_'+targname+'/'+targname+'_allMatchedZPTed.dat',names=True)
    newN = np.genfromtxt(dir+targname+'_absMag.dat',names=True)
    colOs = np.array(oldN.dtype.names)
    colNs = np.array(newN.dtype.names)

    oldCat = np.genfromtxt('catRawMags1305/catDir_'+targname+'/'+targname+'_allMatchedZPTed.dat')
    newCat = np.genfromtxt(dir+targname+'_absMag.dat')

    x_old = np.int(np.where(colOs=='xt1_f606w')[0])
    y_old = np.int(np.where(colOs=='yt1_f606w')[0])

    cs_f606w = np.int(np.where(colOs=='c_star_f606w')[0])
    cs_f814w = np.int(np.where(colOs=='c_star_f814w')[0])

    old_f606w = np.int(np.where(colOs=='magZPT_f606w')[0])
    old_f814w = np.int(np.where(colOs=='magZPT_f814w')[0])

    # x_new = np.int(np.where(colNs=='xt1_f606w')[0])
    # y_new = np.int(np.where(colNs=='yt1_f606w')[0])
    x_new = newN['xt1_f606w']
    y_new = newN['yt1_f606w']

    indOld = np.zeros((len(oldCat),1))
    indOld[:,0] = np.arange(0,len(oldCat),1)

    oldCat = np.hstack((oldCat,indOld))

    master_in = np.zeros((len(x_new),3))
    master_in[:,0] = np.arange(0,len(x_new),1)
    master_in[:,1] = x_new
    master_in[:,2] = y_new

    master_cop = master_in.copy()
    idN_mas,x,y = 0, 1, 2

    cat = oldCat
    idO = len(colOs) # it's in the last column

    matchids = np.zeros((len(master_in),1))

    nF_out = True
    minLen = len(newN)

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,x_old,y_old,idO)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 1
            if matchtol <= 10:
                # master_in = np.hstack((indNew,x_new,y_new))
                master_in = master_cop
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen, 'F606W')

    idN_mas,x_mas,y_mas,idO_mas = 0,1,2,3

    newCols = np.zeros((len(master),4))

    idColNew = master[:,idN_mas]
    idxNew = np.asarray(idColNew,int)
    regNew = newCat[idxNew]

    idColOld = master[:,idO_mas]
    idxOld = np.asarray(idColOld,int)
    regOld = oldCat[idxOld]

    outArr = np.hstack((regNew,regOld[:,[cs_f606w,cs_f814w,old_f606w,old_f814w]]))

    s0 = ' '
    header = s0.join(colNs)
    header += ' c_star_f606w c_star_f814w magZPTold_f606w magZPTold_f814w'

    np.savetxt(dir+'classStarCat_F606Wmatch_'+targname+'_wmag.dat',outArr,header = header)


    return None


def match4cStar814(targname,dir='./',matchtol=3):

    oldN = np.genfromtxt('catRawMags1305/catDir_'+targname+'/'+targname+'_allMatchedZPTed.dat',names=True)
    newN = np.genfromtxt(dir+targname+'_absMag.dat',names=True)
    colOs = np.array(oldN.dtype.names)
    colNs = np.array(newN.dtype.names)

    oldCat = np.genfromtxt('catRawMags1305/catDir_'+targname+'/'+targname+'_allMatchedZPTed.dat')
    newCat = np.genfromtxt(dir+targname+'_absMag.dat')

    x_old = np.int(np.where(colOs=='xt1_f814w')[0])
    y_old = np.int(np.where(colOs=='yt1_f814w')[0])

    cs_f606w = np.int(np.where(colOs=='c_star_f606w')[0])
    cs_f814w = np.int(np.where(colOs=='c_star_f814w')[0])

    # x_new = np.int(np.where(colNs=='xt1_f606w')[0])
    # y_new = np.int(np.where(colNs=='yt1_f606w')[0])
    x_new = newN['xt1_f814w']
    y_new = newN['yt1_f814w']

    indOld = np.zeros((len(oldCat),1))
    indOld[:,0] = np.arange(0,len(oldCat),1)

    oldCat = np.hstack((oldCat,indOld))

    master_in = np.zeros((len(x_new),3))
    master_in[:,0] = np.arange(0,len(x_new),1)
    master_in[:,1] = x_new
    master_in[:,2] = y_new

    master_cop = master_in.copy()

    # master_in = np.hstack((indNew,x_new,y_new))
    idN_mas,x,y = 0, 1, 2

    cat = oldCat
    idO = len(colOs) # it's in the last column

    matchids = np.zeros((len(master_in),1))

    nF_out = True
    minLen = len(newN)

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,x_old,y_old,idO)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 1
            if matchtol <= 10:
                # master_in = np.hstack((indNew,x_new,y_new))
                master_in = master_cop
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen, 'F814W')

    idN_mas,x_mas,y_mas,idO_mas = 0,1,2,3

    newCols = np.zeros((len(master),2))

    idColNew = master[:,idN_mas]
    idxNew = np.asarray(idColNew,int)
    regNew = newCat[idxNew]

    idColOld = master[:,idO_mas]
    idxOld = np.asarray(idColOld,int)
    regOld = oldCat[idxOld]

    outArr = np.hstack((regNew,regOld[:,[cs_f606w,cs_f814w]]))

    s0 = ' '
    header = s0.join(colNs)
    header += ' c_star_f606w c_star_f814w'

    np.savetxt(dir+'classStarCat_F814Wmatch_'+targname+'.dat',outArr,header = header)


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
                master = master[udx]
                matchids_in = matchids_in[udx]
                nF = False

            elif len(udx)==len(master):
                nF = False


    return master,matchids_in
