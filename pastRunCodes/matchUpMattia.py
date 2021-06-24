import numpy as np

psf_dir = '/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf_file = np.genfromtxt(psf_dir + 'HOROLOGIUM_CF.1.TOSEND.CAT',names=True)

file = np.genfromtxt('catRawMags1305/catDir/matchedFLCpsf0906.dat',names=True)
mat_match = np.genfromtxt('pixPosHorIhannah_idmattia.dat',names=True)

mat_g = mat_match[mat_match['nstar']!=-1]
# file_g = file[file['magPSF']>10]
id_col = np.array(mat_g['nstar'],dtype=int)

#  x           y         mag             ra            dec          dx          dy   nstar

# x y mag ra dec dx dy nstar  x y m606c m814x
newArr = np.zeros((len(id_col),12))

for ii, nstar in enumerate(id_col):
    newArr[ii][0] = mat_g['x'][ii]
    newArr[ii][1] = mat_g['y'][ii]
    newArr[ii][2] = mat_g['mag'][ii]
    newArr[ii][3] = mat_g['ra'][ii]
    newArr[ii][4] = mat_g['dec'][ii]
    newArr[ii][5] = mat_g['dx'][ii]
    newArr[ii][6] = mat_g['dy'][ii]
    for jj in range(len(psf_file)):
        if psf_file['nstar'][jj] == nstar:
            newArr[ii][7] = psf_file['x'][jj]
            newArr[ii][8] = psf_file['y'][jj]
            newArr[ii][9] = psf_file['m606c'][jj]
            newArr[ii][10] = psf_file['m814c'][jj]
            newArr[ii][11] = psf_file['nstar'][jj]

header = 'xt1 yt1 mag1 ra dec dx dy xPSF yPSF m606c m814c nstar'
form = '%1.5f %1.5f %1.5f %1.7f %1.7f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %d '

np.savetxt('catRawMags1305/catDir/'+'matchedMattiaInfo.dat',newArr,header=header,fmt=form)
