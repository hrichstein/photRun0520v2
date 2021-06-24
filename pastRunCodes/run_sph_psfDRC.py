import numpy as np
from spherematch import *

dir1='/Volumes/Spare Data/Hannah_Data/mattia/rephotometryquestion/'
psf_file = dir1 + 'HOROLOGIUM_CF.1.TOSEND.CAT'

# dir2 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'
# my_file = dir2+'HOROLOGIUM-I_F606W_at_long.dat'

dir2  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
my_file = dir2 + 'HOROLOGIUM-I_sfErr.dat'

dir3 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'


def run_sph(infile1,infile2,tol=2.78e-4,ra_idx=14,dec_idx=15,suffix='.dat'):

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

    my_header = 'RA_v DEC_v x_v y_v fAper_v fErr_v magAper_v magErr_v magRaw_v magRed_v magAbs_v elong_v ellip_v class_Star_v RA_i DEC_i x_i y_i fAper_i fErr_i magAper_i magErr_i magRaw_i magRed_i magAbs_i elong_i ellip_i class_Star_i corrF_errV corrF_errI corrM_errV corrM_errI'

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

run_sph(psf_file,my_file,suffix='_drc.dat')
matchMags('_drc.dat')
# run_sph(psf_file,my_file,ra_idx=-2,dec_idx=-1,suffix='Dith234.dat')
# matchMags('Dith234.dat')
