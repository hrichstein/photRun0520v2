import numpy as np

file = np.genfromtxt('catRawMags1305/catDir/poorMatch.dat')

RA, DEC,flags,c_star, mag1,mag2,mag3, mag4,xt1, yt1, magZPT,magZPTerr, xPSF_trans, yPSF_trans, xPSF_mas,yPSF_mas,magPSF,id_cat = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ,14 ,15, 16, 17

reg = file[:,[xt1,yt1]]

np.savetxt('catRawMags1305/catDir/poorMatch.reg',reg,fmt='%1.5f %1.5f')
