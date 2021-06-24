import numpy as np
import matplotlib.pyplot as plt


def main():
    targname_arr = np.genfromtxt('./targnamesDirections2.txt',dtype=str)
    se_dir = './catMatchFLCdrc18Oct/seDRCs/'
    for c1,targname in enumerate(targname_arr):
        pu_file = './photUtils20Aug/catDir_' + targname + '/' + targname \
                  + '_filtMatchDRC_pU.dat'
        se_file = se_dir + targname + '_sfErr.dat'
        save_dir = './catMatchFLCdrc18Oct/catDir_' + targname + '/'

        getRefSEpuDRC(targname,seFilename=se_file,puFilename=pu_file,
                      matchtol=3,saveDir=save_dir)


def getRefSEpuDRC(targname,seFilename='None.dat',puFilename='None.dat',
                  matchtol=5,saveDir='./'):

    # The source extractor files will be considered the MASTER here

    # Loading SE files
    seFileN = np.genfromtxt(seFilename, names=True)
    seFile = np.genfromtxt(seFilename)

    # Loading PU Files
    puFileN = np.genfromtxt(puFilename, names=True)
    puFile = np.genfromtxt(puFilename)

    # Getting column names
    seNames = np.array(seFileN.dtype.names)
    puNames = np.array(puFileN.dtype.names)

    # Getting column locations
    xSE = np.int(np.where(seNames=='x_v')[0])
    ySE = np.int(np.where(seNames=='y_v')[0])
    # magSE = np.int(np.where(seNames=='magRaw_v')[0])

    xPU = np.int(np.where(puNames=='xcenter_f606w')[0])
    yPU = np.int(np.where(puNames=='ycenter_f606w')[0])
    # magPU = np.int(np.where(puNames=='magr_f606w')[0])

    # Making an id column for the matching array
    idCol = np.zeros((len(puFile),1))
    idCol[:,0] = np.arange(0,len(puFile),1)

    puFile_wid = np.hstack((puFile,idCol))
    id_idx = len(puNames)

    # Getting the 50 brightest stars
    s50 = np.argsort(seFileN['magRaw_v'])[:50]
    se_50 = seFile[s50]

    p50 = np.argsort(puFileN['magr_f606w'])[:50]
    pu_50 = puFile_wid[p50]

    master_in = se_50[:,[xSE,ySE]]
    xse_idx, yse_idx = 0, 1

    cat = pu_50

    matchtol = matchtol
    nF_out = True

    while nF_out:
        master, matchids = matchlistID(master_in,cat,matchtol,xse_idx,yse_idx,
                                       xPU,yPU,id_idx)

        if len(master)>=int(6):  # Need 6 stars to do the 6D transform
            nF_out = False
            print('Minimum Number Reached:{0:d}'.format(len(master)),targname)

        else:
            print('Need more stars...')
            master_in = se_50[:,[xSE,ySE]]  # need to reset master_in just
            # in case it was changed from the function
            matchtol += 1  # Increasing the match tolerance.
            # Change depending on the step size you need to take

    # Combining the master x,y with the matching indices for the cat array
    master = np.hstack((master,matchids))
    xse_idx, yse_idx, id_match = 0, 1, 2

    # Using the indices (have to do a little type-magic first)
    idxCol = master[:,id_match]
    idxP = np.asarray(idxCol,int)
    regP = puFile_wid[idxP]

    # Since these are just reference stars, only need the x,y info
    newCols = np.zeros((len(master),2))
    newCols[:,0] = regP[:,xPU]
    newCols[:,1] = regP[:,yPU]

    outArr = np.hstack((master,newCols))

    header = 'xSE ySE idPU xPU yPU'
    form = '%1.5f %1.5f %d %1.5f %1.5f'

    outName = saveDir + 'sePUdrcRef'
    np.savetxt(outName + '.dat', outArr, header=header, fmt=form)

    # Prep for plotting, then plotting
    xse, yse, idP, xpu, ypu = 0, 1, 2, 3, 4

    fig, ax = plt.subplots(figsize=(4,4))

    ax.scatter(outArr[:,xse],outArr[:,yse],label='SE Pos',s=50,color='black')
    ax.scatter(outArr[:,xpu],outArr[:,ypu],label='PU Pos',s=10,color='magenta')

    ax.legend()
    ax.set_title(targname)

    plt.savefig(outName + '.png',dpi=600,bbox_inches='tight')
    plt.close()  # !! Very important line!!!

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
            # udx is the array of unique indices. I see if this is less than
            # the master length. If so, I use udx to get the relevant, unique
            # sources.
            u, udx = np.unique(matchids_in,return_index=True)

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


if __name__ == '__main__':
    main()
#
