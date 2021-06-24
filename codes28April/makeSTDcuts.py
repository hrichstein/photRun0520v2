import numpy as np


def makeSTDcuts(filt,dir='./'):

    dataN = np.genfromtxt(dir+'matched_flcs_'+filt+'.dat',names=True)
    data = np.genfromtxt(dir+'matched_flcs_'+filt+'.dat')

    colNs = np.array(dataN.dtype.names)

    mag1 = np.int(np.where(colNs=='magr1')[0])
    mag2 = np.int(np.where(colNs=='magr2')[0])
    mag3 = np.int(np.where(colNs=='magr3')[0])
    mag4 = np.int(np.where(colNs=='magr4')[0])

    newCol = np.zeros((len(data),2))

    for dd in range(len(data)):
        idx = np.array([mag1,mag2,mag3,mag4],dtype=int)  # columns of mags
        std = np.nanstd(data[dd][idx])

        newCol[dd,0] = np.nanmean(data[dd][idx])
        newCol[dd,1] = std/np.sqrt(len(idx))

    data = np.hstack((data, newCol))

    s0 = ' '
    header = s0.join(colNs)

    header += ' mean err'

    np.savetxt(dir+'magSTDcut_'+filt+'.dat',data,header=header)

    # print("Percentage Kept",kept/total)

    return None
