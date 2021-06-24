import numpy as np
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord

file = np.genfromtxt('../Hannah_Data/drcTargInfo_new.dat',names=True)
colNs = np.array(file.dtype.names)

fileCat = np.genfromtxt('../Hannah_Data/drcTargInfo_new.dat',dtype=str)

targ = np.int(np.where(colNs=='TARGNAME')[0])
filt1 = np.int(np.where(colNs=='FILTER1')[0])
filt2 = np.int(np.where(colNs=='FILTER2')[0])

ra = file['RA']
dec = file['DEC']

v_ebv=[[] for ll in range(len(file))]
i_ebv=[[] for ll in range(len(file))]

for tt in range(len(file)):
    coords = SkyCoord(ra[tt],dec[tt],frame='icrs',unit='degree')
    sfd = SFDQuery()
    v_ebv[tt] = 2.488 * sfd(coords)
    i_ebv[tt] = 1.536 * sfd(coords)

v_ebv = np.array(v_ebv)
i_ebv = np.array(i_ebv)

outArr = np.zeros((len(file),4))

# outArr[:,0] = fileCat[:,targ]
# outArr[:,1] = fileCat[:,filt1]
# outArr[:,2] = fileCat[:,filt2]
outArr[:,0] = v_ebv
outArr[:,1] = i_ebv
outArr[:,2] = ra
outArr[:,3] = dec

out = np.hstack((fileCat[:,[targ,filt1,filt2]],outArr))

header = 'TARGNAME FILTER1 FILTER2 V_EBV I_EBV RA DEC'

np.savetxt('sfdCorrection.dat',out,header=header,fmt='%s')
