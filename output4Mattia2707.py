import numpy as np

file = np.genfromtxt('catRawMags1305/catDir_SAGITTARIUS-II/sgrExamplePts.dat',names=True)

outArr = np.zeros((len(file),6))

outArr[:,0] = file['xR']
outArr[:,1] = file['yR']
outArr[:,2] = file['magrF']
outArr[:,3] = file['magDcorrNZPT']
outArr[:,4] = file['dCorr']
outArr[:,5] = file['magrF'] - 26.779 - (2.5*np.log10(1138.)) - (2.5*np.log10(0.770))

header= 'x_raw y_raw mag_raw mag_corrected mag_correction mag_instr_raws'

np.savetxt('catRawMags1305/catDir_SAGITTARIUS-II/sgr2ex_Mattia.dat',outArr,header=header)
