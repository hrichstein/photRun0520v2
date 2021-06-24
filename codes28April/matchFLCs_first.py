import numpy as np


def matchFLCs(targname,filt,jdanUse,dir='./',matchtol=2.5):

    masterN = np.genfromtxt(dir+jdanUse[0]+'_'+filt+"_pu.dat",names=True)
    master = np.genfromtxt(dir+jdanUse[0]+'_'+filt+"_pu.dat")

    colNs = np.array(masterN.dtype.names)

    xt = np.int(np.where(colNs=='xDRC')[0])
    yt = np.int(np.where(colNs=='yDRC')[0])
    id = np.int(np.where(colNs=='id')[0])

    matchids = np.zeros((len(master), (len(jdanUse)-1)),dtype=int)

    for dd in range(len(jdanUse)-1):
        # Load catalogs
        cat = np.genfromtxt(dir+jdanUse[dd+1]+'_'+filt+"_pu.dat")

        master, matchids = matchlistID(master,cat,matchids,xt,yt,xt,yt,id,dd,
                                       matchtol)

    jj = 0
    while jj < len(jdanUse):
        cat = np.genfromtxt(dir+jdanUse[jj]+'_'+filt+"_pu.dat")

        if jj==0:
            idcol = master[:,id]
        else:
            idcol = matchids[:,jj-1]

        idx = np.asarray(idcol,int)
        reg = cat[idx]

        if jj==0:
            tempArr = reg
        else:
            tempArr = np.hstack((tempArr,reg))

        jj += 1

    outArr = tempArr

    header1 = '1 '.join(colNs)
    header2 = '2 '.join(colNs)
    header3 = '3 '.join(colNs)
    header4 = '4 '.join(colNs)

    header = header1 + '1 ' + header2 + '2 ' + header3 + '3 ' + header4 + '4'

    np.savetxt(dir + 'matched_flcs_' + filt + '.dat',outArr,header=header)

    return None


def matchlistID(master,cat,matchids_in,x1,y1,x2,y2,id_mat,dd=0,matchtol=2.5):

    nF = True
    row = 0

    while nF:

        matchrows = cat[(abs(master[row][x1] - cat[:,x2])
                        <= matchtol) & (abs(master[row][y1] - cat[:,y2])
                        <= matchtol)]

        if (len(matchrows) == 1):
            matchids_in[row][dd] = matchrows[0][id_mat]
            row += 1

        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for ee in range(len(matchrows)):
                distDiff[ee] = np.sqrt((master[row][x1]
                                       - matchrows[ee][x2])**2
                                       + (master[row][y1]
                                       - matchrows[ee][y2])**2)
            small = np.argmin(distDiff)
            matchids_in[row][dd] = matchrows[small][id_mat]
            row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            nF = False

    return master,matchids_in

#
