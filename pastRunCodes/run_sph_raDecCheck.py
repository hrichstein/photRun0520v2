import numpy as np
from spherematch import *

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf_file = dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT'

dir2 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'
my_file = dir2+'HOROLOGIUM-I_F606W_at_long.dat'

print(len(psf_file))
print(len(my_file))

def run_sph(infile1,infile2,tol=2.78e-4,ra_idx=0,dec_idx=1,suffix='Dith1.dat'):

    x1 = np.loadtxt(psf_file, dtype='float')
    x2 = np.loadtxt(my_file, dtype='float')

    ra1 = x1[:, -2]
    dec1 = x1[:, -1]
    ra2 = x2[:, ra_idx]
    dec2 = x2[:, dec_idx]

    idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tol=tol, nnearest=1)

    outfile1 = dir2+'psf_match_'+suffix
    np.savetxt(outfile1, idx1, fmt='%4i')
    outfile2 = dir2+'my_match_'+suffix
    np.savetxt(outfile2, idx2, fmt='%4i')

    print(len(x1))
    print(len(x2))

    return None

def matchMags(suffix='Dith1.dat'):

    my_match = dir2+'my_match_'+suffix
    psf_match = dir2+'psf_match_'+suffix

    filters = ['F606W','F814W']

    my_header = 'wcsRA wcsDEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1 xr2 yr2 xc2 yc2 xo2 yo2 xr3 yr3 xc3 yc3 xo3 yo3 xr4 yr4 xc4 yc4 xo4 yo4 xt2 yt2 wra2 wdec2 xt3 yt3 wra3 wdec3 xt4 yt4 wra4 wdec4 meanRA_234 meanDEC_234'

    # my_form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f'

    psf_header = "x y m606c m814c nstar sat606 sat814 camera m606 s606 q606 o606 f606 g606 rxs606 sky606 rmssky606 m814 s814 q814 o814 f814 g814 rxs814 sky814 rmssky814 ra dec"

    for ff in range(len(filters)):
        dataFile = np.genfromtxt(my_file)
        idxFile = np.genfromtxt(my_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir2+'my'+suffix,outArr,header=my_header)

        dataFile = np.genfromtxt(psf_file)
        idxFile = np.genfromtxt(psf_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir2+'psf'+suffix,outArr,header=psf_header)

    return None

run_sph(psf_file,my_file,suffix='Dith1.dat')
matchMags('Dith1.dat')
run_sph(psf_file,my_file,ra_idx=-2,dec_idx=-1,suffix='Dith234.dat')
matchMags('Dith234.dat')
