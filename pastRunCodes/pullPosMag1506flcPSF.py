import numpy as np

dir = 'catRawMags1305/catDir/'

def pullFunc(filter):
    if filter=='F606W':
        master= np.genfromtxt(dir+'flcPSF_idx_1506_606.dat')
        cat = np.genfromtxt(dir+'flcPSFpos1106.dat')

        RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19, 20, 21

        outname= 'matchedFLCpsf1506_606wMean.dat'

    elif filter=='F814W':
        master= np.genfromtxt(dir+'flcPSF_idx_1506_814.dat')
        cat = np.genfromtxt(dir+'flcPSFpos0906.dat')

        RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xt1, yt1, mean, stdev, cut_flag, idx_cut, num_abv_std, magZPT, magZPTerr,xPSF_trans, yPSF_trans, id_new = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19

        outname= 'matchedFLCpsf1506_814wMean.dat'

    x,y,magr,magerr,nstar,id_cat = 0,1,2,3,4,5

    newCols = np.zeros((len(master),9))

    idCol = master[:,id_cat]
    idx = np.asarray(idCol,int)

    reg = cat[idx]

    newCols[:,0] = reg[:,RA]
    newCols[:,1] = reg[:,DEC]
    newCols[:,2] = reg[:,flags]
    newCols[:,3] = reg[:,c_star]
    newCols[:,4] = reg[:,magZPT]
    newCols[:,5] = reg[:,magZPTerr]
    newCols[:,6] = reg[:,xPSF_trans]
    newCols[:,7] = reg[:,yPSF_trans]
    newCols[:,8] = reg[:,mean]


    outArr = np.hstack((newCols,master[:,[x,y,magr,magerr,nstar,id_cat]]))
    header= 'RA DEC flags c_star magZPT magZPTerr xPSF_trans yPSF_trans mean xPSF_mas yPSF_mas magPSF magErrPSF nstar id_cat'

    np.savetxt(dir+outname,outArr,header=header)

    return None

pullFunc('F606W')
pullFunc('F814W')
