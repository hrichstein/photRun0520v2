import numpy as np

dir2 = 'catRawMags1305/catDir/'

xPSF, yPSF, m606cPSF, m814cPSF, s606PSF, s814PSF, nstarPSF, nstarAPER, idAPER, xAPER, yAPER, m606cAPER, m814cAPER, s606APER, s814APER, sky606APER, sky814APER, ra, dec, id =  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19

aperCat = np.genfromtxt(dir2+'matchedPSFaper2906_tc.dat')

sort_idx = np.argsort(aperCat[:,m606cPSF])[:50]

idxCol = np.zeros((len(sort_idx),1))
idxCol[:,0] = np.arange(0,len(sort_idx),1)

psf = aperCat[sort_idx]
psf = np.hstack((psf,idxCol))

# np.savetxt('psf50_bkgd.reg',psf[:,[ra,dec]],fmt='%1.7f')

np.savetxt('psf50_bkgd.dat',psf[:,[ra,dec,xPSF,yPSF,id]],fmt='%1.7f %1.7f %1.5f %1.5f %d')

sort_idx = np.argsort(aperCat[:,m606cAPER])[:50]

aper = aperCat[sort_idx]
aper = np.hstack((aper,idxCol))

# np.savetxt('aper50_bkgd.reg',aper[:,[ra,dec]],fmt='%1.7f')

np.savetxt('aper50_bkgd.dat',aper[:,[ra,dec,xAPER,yAPER,id]],fmt='%1.7f %1.7f %1.5f %1.5f %d')
