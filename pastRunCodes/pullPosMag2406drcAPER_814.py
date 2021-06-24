import numpy as np

dir = 'catRawMags1305/catDir/'

master= np.genfromtxt(dir+'drcAPER_idx_2406_f814w.dat')
cat = np.genfromtxt(dir+'drc2APER_f814w_t.dat')

RA_v, DEC_v, x_v, y_v, fAper_v, fErr_v, magAper_v, magErr_v, magRaw_v, magRed_v, magAbs_v, elong_v, ellip_v, class_Star_v, RA_i, DEC_i, x_i, y_i, fAper_i, fErr_i, magAper_i, magErr_i, magRaw_i, magRed_i, magAbs_i, elong_i, ellip_i, class_Star_i, corrF_errV, corrF_errI, corrM_errV, corrM_errI,xt_i,yt_i,id = 0, 1, 2 ,3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34

x,y,magr,id_cat = 0,1,2,3

newCols = np.zeros((len(master),4))

idCol = master[:,id_cat]
idx = np.asarray(idCol,int)

reg = cat[idx]

newCols[:,0] = reg[:,magRaw_i]
newCols[:,1] = reg[:,corrM_errI]
newCols[:,2] = reg[:,xt_i]
newCols[:,3] = reg[:,yt_i]


outArr = np.hstack((newCols,master[:,[x,y,magr,id_cat]]))
header= 'mag814 magErr814 xt_i yt_i xAPER yAPER magAPER idAPER'

np.savetxt(dir+'matchedDRCaper2406_814.dat',outArr,header=header)
