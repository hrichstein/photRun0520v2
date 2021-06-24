import numpy as np

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
