import numpy as np
from spherematch import *

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf_file = dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT'

# dir2 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'
# my_file = dir2+'HOROLOGIUM-I_F606W_at_long.dat'

dir2  = '/Volumes/Spare Data/Hannah_Data/hor1dir2804/'
my_file = dir2 + 'HORI_comb0605.dat'

dir3 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'


def run_sph(infile1,infile2,tol=2.78e-4,ra_idx=11,dec_idx=12,suffix='.dat'):

    x1 = np.loadtxt(psf_file, dtype='float')
    x2 = np.loadtxt(my_file, dtype='float')

    ra1 = x1[:, -2]
    dec1 = x1[:, -1]
    ra2 = x2[:, ra_idx]
    dec2 = x2[:, dec_idx]

    idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tol=tol, nnearest=1)

    outfile1 = dir3+'psf_match_'+suffix
    np.savetxt(outfile1, idx1, fmt='%4i')
    outfile2 = dir3+'my_match_'+suffix
    np.savetxt(outfile2, idx2, fmt='%4i')

    print(len(x1))
    print(len(x2))

    return None

def matchMags(suffix='.dat'):

    my_match = dir3+'my_match_'+suffix
    psf_match = dir3+'psf_match_'+suffix

    filters = ['F606W','F814W']

    my_header = 'RA_f606w DEC_f606w flux_f606w c_star_f606w x_f606w y_f606w rawMean_f606w zptMean_f606w rawErr_f606w zptErr_f606w  pos_std_f606w RA_f814w DEC_f814w flux_f814w c_star_f814w x_f814w y_f814w rawMean_f814w zptMean_f814w rawErr_f814w zptErr_f814w  pos_std_f814w'

    # my_form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f'

    psf_header = "x y m606c m814c nstar sat606 sat814 camera m606 s606 q606 o606 f606 g606 rxs606 sky606 rmssky606 m814 s814 q814 o814 f814 g814 rxs814 sky814 rmssky814 ra dec"

    for ff in range(len(filters)):
        dataFile = np.genfromtxt(my_file)
        idxFile = np.genfromtxt(my_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir3+'my'+suffix,outArr,header=my_header)

        dataFile = np.genfromtxt(psf_file)
        idxFile = np.genfromtxt(psf_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir3+'psf'+suffix,outArr,header=psf_header)

    return None

run_sph(psf_file,my_file,suffix='.dat')
matchMags('.dat')
# run_sph(psf_file,my_file,ra_idx=-2,dec_idx=-1,suffix='Dith234.dat')
# matchMags('Dith234.dat')
