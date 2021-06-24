import numpy as np

dir = 'catRawMags1305/catDir/'

# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'psfAPER_idx_1706.dat')

psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'

psf = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT')
aper = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.2.TOSEND.CAT')

x, y, m606c, m184c, s606, s814, nstarPSF, nstarAPER, idAPER = 0, 1, 2, 3, 4, 5, 6, 7, 8

newCols = np.zeros((len(master),10))

idCol = master[:,idAPER]
idx = np.asarray(idCol,int)
reg = aper[idx]

x, y, m606c, m814c, nstar, sat606, sat814, camera, m606, s606, q606, o606, f606, g606, rxs606, sky606, rmssky606, m814, s814, q814, o814, f814, g814, rxs814, sky814, rmssky814, ra, dec = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 #psf

newCols[:,0] = reg[:,x]
newCols[:,1] = reg[:,y]
newCols[:,2] = reg[:,m606c]
newCols[:,3] = reg[:,m814c]
newCols[:,4] = reg[:,s606]
newCols[:,5] = reg[:,s814]
newCols[:,6] = reg[:,sky606]
newCols[:,7] = reg[:,sky814]
newCols[:,8] = reg[:,ra]
newCols[:,9] = reg[:,dec]

##################################

outArr = np.hstack((master,newCols))
header= 'xPSF yPSF m606cPSF m814cPSF s606PSF s814PSF nstarPSF nstarAPER idAPER xAPER yAPER m606cAPER m814cAPER s606APER s814APER sky606APER sky814APER ra dec'

# np.savetxt(dir+'matchedPSFpsf1106.dat',outArr,header=header)
np.savetxt(dir+'matchedPSFaper2906.dat',outArr,header=header)
