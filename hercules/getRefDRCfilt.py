"""
Getting the reference stars between the two filters so that I can do a 6D
linear transform.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# IMPORTANT VARIABLES FOR THIS PROCESS
work_dir = './'  # this goes to photRun0520
# targFile = work_dir + 'targnamesDirections2.txt'
dateDef = '18Mar'
suffix_ = '_pu.dat'


def getRefDRCfilt(targname,f6d,f8d,dir='./',matchtol=3,suffix='_pu.dat',
                  setnum=0):
    """
    Replace F606W and F814W (and 606 and 814) throughout this if using
    different filters
    """

    x = os.listdir(dir)
    for ii in x:
        if ii.endswith(f6d+'_drc_F606W' + suffix):
            f606w_file = ii
        elif ii.endswith(f8d+'_drc_F814W' + suffix):
            f814w_file = ii

    f606wN = np.genfromtxt(dir+f606w_file,names=True)
    f606w = np.genfromtxt(dir+f606w_file)

    f814wN = np.genfromtxt(dir+f814w_file,names=True)
    f814w = np.genfromtxt(dir+f814w_file)

    col606 = np.array(f606wN.dtype.names)
    # col814 = np.array(f814wN.dtype.names)  # names will be the same between
    # filters

    # Getting indices of columns

    # Will be the same for both filters
    x = np.int(np.where(col606=='xcenter')[0])
    y = np.int(np.where(col606=='ycenter')[0])
    mag = np.int(np.where(col606=='magr')[0])

    idCol = np.int(np.where(col606=='id')[0])  # the row index

    # Sorting to get the 50 brightest stars
    v50 = np.argsort(f606wN['magr'])[:50]  # using v abbreviation because
    # F606W is a broad V-band
    f606w_50 = f606w[v50]

    i50 = np.argsort(f814wN['magr'])[:50]  # using i abbrev. b/c F814W is
    # broad I-band
    f814w_50 = f814w[i50]

    master_in = f606w_50[:,[idCol,x,y,mag]]  # selecting the columns I need
    idV,xv,yv,magv = 0,1,2,3  # associating short vars with indices

    cat = f814w_50  # this has ALL of the columns from F814W catalog
    # matchids = np.zeros((len(master_in),1))

    nF_out = True
    # matchtol=matchtol

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,xv,yv,x,y,idCol)

        if len(master)>=int(6):  # because it's a 6D transformation
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)
        else:
            print('Need More Stars')
            master_in = f606w_50[:,[idCol,x,y,mag]]  # resetting, just in case
            # matchids = np.zeros((len(master_in),1))
            matchtol += 1

    master = np.hstack((master,matchids))

    idV, xv, yv, magv, idI = 0, 1, 2, 3, 4  # index columns of master

    newCols = np.zeros((len(master),3))

    idxCol = master[:,idI]
    idxI = np.asarray(idxCol,int)
    regI = f814w[idxI]

    newCols[:,0] = regI[:,x]
    newCols[:,1] = regI[:,y]
    newCols[:,2] = regI[:,mag]

    outArr = np.hstack((master,newCols))

    idV, xv, yv, magv, idI, xi, yi, magi = 0, 1, 2, 3, 4, 5, 6, 7

    header = 'id_f606w x_f606w y_f606w magr_f606w id_f814w x_f814w y_f814w'
    header += ' magr_f814w'  # had to take add line since first was too long
    form = '%d %1.6f %1.6f %1.6f %d %1.6f %1.6f %1.6f'

    outName = dir+'drcFiltRef_'+targname+str(setnum)
    # np.savetxt(outName+'.dat',outArr,header=header,fmt=form)
    np.savetxt(outName+'.dat',outArr,header=header,fmt=form)

    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(outArr[:,xv],outArr[:,yv],label='F606W',s=50,color='black')
    ax.scatter(outArr[:,xi],outArr[:,yi],label='F814W',s=20,color='magenta')

    ax.legend()
    ax.set_title(targname)

    # plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')
    plt.savefig(outName+'.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None


def matchlistID(master,cat,matchtol,x1,y1,x2,y2,id_mat):

    """
    Input:
    master: an array, n by n, with the relevant information in columns
    cat: the array that will be searched for matches, n by n with relevant
    info in columns
    x1: the column index where the x-coordinates are listed in the master array
    y1: the column index where the y-coordinates are listed in the master array
    x2: the column index where the x-coordinates are listed in the "cat"
    (matching) array
    y2: the column index where the y-coordinates are listed in the "cat"
    (matching) array
    id_mat: column in the matching array with the indices of sources as listed
    in the raw file

    Output:
    master: a shortened version of the master array that was input;
    only sources in the master than have a match in the cat are listed
    match_ids: an len(master) by 1 array that has the indices of the best-
    matching sources from the input cat array (ids come from the id column of
    the cat array; if the cat array was shortened from a longer cat array (to
    which this index list will be applied, the indices should have come from
    that longer cat array)). I'm matching based on index position, although in
    the future, to make things easier (possibly?) use real (kind of random,
    not necessarily listed order) IDs and match based on what the column ID
    says
    """

    matchids_in = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:
        # Finding the difference in the x and y positions between the master
        # array, row by row
        # Checking that the difference in x and y (separately) is smaller than
        # the match tolerance
        # The statements in the brackets return indices for which the
        # conditions are true
        # Then, matchrows is a listing of the rows in the cat array where the
        # sources meet the positional requirements.
        matchrows = cat[(abs(master[row][x1] - cat[:,x2])
                        <= matchtol) & (abs(master[row][y1] - cat[:,y2])
                        <= matchtol)]

        # If only one source met the tolerance criteria, the index value for
        # that source is put into the matchids array. It will go in the row
        # corresponding to the row where the master source was.
        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[0][id_mat]
            row += 1

        # If there is more than one source that meets the criteria,
        # I calculate the distance between all of the match sources
        # and the master source. I put these in an array and find the (row)
        # index of the minimum distance. I proceed to put the cat array
        # source index into the matchids array.
        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt((master[row][x1]
                                       - matchrows[dd][x2])**2
                                       + (master[row][y1]
                                       - matchrows[dd][y2])**2)
            small = np.argmin(distDiff)
            matchids_in[row][0] = matchrows[small][id_mat]
            row += 1

        # If there is nothing that meets the criteria, the master source
        # row is removed, as well as the corresponding row in the matchids
        # array.
        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        # If the row counter is longer than the length of the master,
        # we've reached the end of the distance tabulations. I do a uniqueness
        # check as if there's a repeat in the match_ids, it means multiple
        # master sources matched with the same cat source. (I'm going
        # methodically through the master list, but not removing best-matching
        # sources from the cat array.)
        if (row >= len(master)):
            u, udx = np.unique(matchids_in,return_index=True)
            # udx is the array of unique indices. I see if this is less than
            # the master length. If so, I use udx to get the relevant, unique
            # sources.
            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]

                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(
                    matchtol,len(master)))
                nF = False

            elif len(udx)==len(master):
                print("Pixel Tolerance: {0:d}, Number Stars: {1:d}".format(
                    matchtol,len(master)))
                nF = False
        # A different way to deal with this would be to invert the matching
        # algorithm. Take the two matched lists and start running through the
        # cat array side. The cat sources with multiple master sources should
        # be fixed by taking the master source with the smallest distance. As
        # it is now, I just remove all the sources associated with repeats.

    return master,matchids_in


#
