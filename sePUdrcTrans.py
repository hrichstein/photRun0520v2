""" This transforms PU positions to SE """

import numpy as np
import matplotlib.pyplot as plt

from linear6d import test_linear


def main():
    targname_arr = np.genfromtxt('./targnamesDirections2.txt', dtype=str)
    for c1, targname in enumerate(targname_arr):
        ref_file = './catMatchFLCdrc18Oct/catDir_'+ targname \
                   + '/sePUdrcRef.dat'
        all_file = './photUtils20Aug/catDir_'+ targname + '/' + targname \
                   + '_filtMatchDRC_pU.dat'
        save_dir = './catMatchFLCdrc18Oct/catDir_'+ targname + '/'

        linTransSEpu(targname, refFile=ref_file, allFile=all_file,
                     saveDir=save_dir)


def linTransSEpu(targname, refFile='None.dat', allFile='None.dat',
                 saveDir='./'):

    # Going from photutils to source extractor; although it shouldn't
    # really matter which direction I go as long as the eventual
    # matching code knows what I'm doing here

    # Loading the reference file and putting the names in an array;
    # the "col" randomly refers to column.
    # The -N suffix means that array was loaded with Names,
    # and is therefore not has a form (n,)
    mFileN = np.genfromtxt(refFile, names=True)
    mFile = np.genfromtxt(refFile)
    colMs = np.array(mFileN.dtype.names)

    # Loading the file with source positions to be transformed;
    # calling this the all file because it has all the points,
    # not just the few reference points
    aFileN = np.genfromtxt(allFile, names=True)
    aFile = np.genfromtxt(allFile)
    colAs = np.array(aFileN.dtype.names)

    # Putting the photutils positions into the match array
    # to be used as references
    match_arr = np.zeros((len(mFile), 2))
    xtrans = np.int(np.where(colMs=='xPU')[0])
    ytrans = np.int(np.where(colMs=='yPU')[0])

    match_arr[:,0] = mFile[:,xtrans]
    match_arr[:,1] = mFile[:,ytrans]

    # Putting the source extractor positions into the master array
    master_arr = np.zeros((len(mFile),2))
    x1 = np.int(np.where(colMs=='xSE')[0])
    y1 = np.int(np.where(colMs=='ySE')[0])

    master_arr[:,0] = mFile[:,x1]
    master_arr[:,1] = mFile[:,y1]

    # Weights are needed for the 6D code, but we've
    # never really taken them into account
    weights = np.zeros((len(master_arr)))
    weights.fill(1.0)

    # Filling the array that will have the values that need
    # to be transformed; in this case, they are going from PU to SE;
    # variable names with _bt refer to "before tranform"
    all_arr = np.zeros((len(aFile), 2))
    x_bt = np.int(np.where(colAs=='xcenter_f606w')[0])
    y_bt = np.int(np.where(colAs=='ycenter_f606w')[0])

    all_arr[:,0] = aFile[:,x_bt]
    all_arr[:,1] = aFile[:,y_bt]

    s0 = ' '
    header = s0.join(colAs)  # "pulling" the header from the original file
    header += ' xSEf606w ySEf606w'  # adding new column names; make sure
    # there's a SPACE at the beginning, so the names don't run together

    outName = saveDir + 'puWtrans'

    # Call to the 6D code. More info in that script on the inputs;
    # the above portion of code put together the necessary arrays;
    # new_match is the positions of the original reference stars
    # transformed into the new frame
    # new_all is the list of all the transformed positions
    new_match, new_all = test_linear(match_arr[:,0], match_arr[:,1],
                                     master_arr[:,0], master_arr[:,1],
                                     weights, weights, all_arr[:,0],
                                     all_arr[:,1])

    # Tacking the columns onto the original array
    outArr = np.hstack((aFile, new_all))

    np.savetxt(outName + '.dat', outArr, header=header)

    # Making a plot of the transformed reference stars to see
    # if anything went wrong
    makePlot(targname,match_arr[:,0],match_arr[:,1],
             new_match[:,0],new_match[:,1],master_arr[:,0], master_arr[:,1],
             label_1='Original in PU',label_2='New in PU 2 SE',
             label_3='Original in SE',outname=outName)

    return None


def makePlot(targname,x1,y1,x2,y2,x3,y3,label_1,
             label_2,label_3,outname):

    fig, ax = plt.subplots(figsize=(4,4))

    ax.scatter(x3,y3,label=label_3,s=50,color='black')
    ax.scatter(x1,y1,label=label_1,s=20,color='magenta')
    ax.scatter(x2,y2,label=label_2,s=10,color='dodgerblue')

    ax.legend()
    ax.set_title(targname)

    plt.savefig(outname + '_matchCheck.png',dpi=600,bbox_inches='tight')
    plt.close()

    return None


if __name__ == '__main__':
    main()

#
