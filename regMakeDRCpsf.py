import numpy as np

workDir='./catRawMags1305/catDir/'
fileDir = './catRawMags1305/'


# file = 'jdan21l8q_HOROLOGIUM-I_F814W_wMag.dat'

# data = np.genfromtxt(fileDir+file)

# flags, RA, DEC, xr, yr, flux, c_star, magr, id = 0, 1, 2, 3, 4, 5, 6, 7, 8


# dith1idx = np.argsort(data[:,magr])[:1000]
# dith1_1000 = data[dith1idx]

# np.savetxt(workDir+'HOROLOGIUM-I_F814W_dith1_1000.reg',dith1_1000[:,[RA,DEC]],fmt='%1.6f')
# dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
# psf_file = np.genfromtxt(dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT')
#
dir3  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
drc_file = np.genfromtxt(dir3 + 'HOROLOGIUM-I_sfErr.dat')

RA_v, DEC_v, x_v, y_v, fAper_v, fErr_v, magAper_v, magErr_v, magRaw_v, magRed_v, magAbs_v, elong_v, ellip_v, class_Star_v, RA_i, DEC_i, x_i, y_i, fAper_i, fErr_i, magAper_i, magErr_i, magRaw_i, magRed_i, magAbs_i, elong_i, ellip_i, class_Star_i, corrF_errV, corrF_errI, corrM_errV, corrM_errI,id = 0, 1, 2 ,3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
#
# newCol = np.zeros((len(psf_file),1),dtype=int)
# newCol[:,0] = np.arange(0,len(psf_file),1)
# id = int(32)
drc_idx = np.argsort(drc_file[:,magRaw_i])[:100]
#
drc_1000 = drc_file[drc_idx]
#
newCol = np.zeros((len(drc_1000),1),dtype=int)
newCol[:,0] = np.arange(0,len(drc_1000),1)

drc_1000 = np.hstack((drc_1000,newCol))
# print(len(drc_file))
# print(len(drc_1000))
#
np.savetxt(workDir+'HOROLOGIUM-I_F814W_drc100_1706.reg',drc_1000[:,[RA_i,DEC_i]],fmt='%1.6f')

np.savetxt(workDir+'HOROLOGIUM-I_F814W_drc100_1706.dat',drc_1000[:,[RA_i,DEC_i,x_i,y_i,id]],fmt='%1.6f %1.6f %1.5f %1.5f %d')

# x, y, m606c, m814c, nstar, sat606, sat814, camera, m606, s606, q606, o606, f606, g606, rxs606, sky606, rmssky606, m814, s814, q814, o814, f814, g814, rxs814, sky814, rmssky814, ra, dec,id = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,28
#
# newCol = np.zeros((len(psf_file),1),dtype=int)
# newCol[:,0] = np.arange(0,len(psf_file),1)
#
# psf_file = np.hstack((psf_file,newCol))
#
# psf_g = psf_file[psf_file[:,m814c]>10]
#
# psf_idx = np.argsort(psf_g[:,m814c])[:100]
# psf_100 = psf_g[psf_idx]
# #
# # # psf_out = np.array([psf_50[10],psf_50[12],psf_50[13],psf_50[15],psf_50[20],psf_50[25]])
# #
# # header = 'x y m606c m814c nstar sat606 sat814 camera m606 s606     q606 o606 f606 g606 rxs606 sky606 rmssky606 m814 s814 q814 o814 f814 g814 rxs814 sky814 rmssky814 ra dec id'
# #
# # np.savetxt(workDir+'HOROLOGIUM-I_F814W_psf50.reg',psf_100[:,[ra,dec]],fmt='%1.6f')
# header= 'x y ra dec id'
#
# np.savetxt(workDir+'HOROLOGIUM-I_F606W_psf100ID.dat',psf_100[:,[x,y,ra,dec,id]],header=header,fmt='%1.5f %1.5f %1.7f %1.7f %d')
