import numpy as np
from spherematch import *

dir1  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
drc_file = dir1 + 'HOROLOGIUM-I_sfErr.dat'

# dir2 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'
# my_file = dir2+'HOROLOGIUM-I_F606W_at_long.dat'

dir2  = '/Volumes/Spare Data/Hannah_Data/hor1dir2804/'
my_file = dir2 + 'HORI_comb0605.dat'

dir3 = '/Volumes/Spare Data/photRun0520/catRawMags1305/catDir/'


def run_sph(infile1,infile2,tol=2.78e-4,ra_idx=14,dec_idx=15,suffix='.dat'):

    x1 = np.loadtxt(drc_file, dtype='float')
    x2 = np.loadtxt(my_file, dtype='float')

    ra1 = x1[:, ra_idx]
    dec1 = x1[:, dec_idx]
    ra2 = x2[:,11]
    dec2 = x2[:, 12]

    idx1, idx2, ds = spherematch(ra1, dec1, ra2, dec2, tol=tol, nnearest=1)

    outfile1 = dir3+'drc_match_'+suffix
    np.savetxt(outfile1, idx1, fmt='%4i')
    outfile2 = dir3+'my_match_'+suffix
    np.savetxt(outfile2, idx2, fmt='%4i')

    print(len(x1))
    print(len(x2))

    return None

def matchMags(suffix='.dat'):

    my_match = dir3+'my_match_'+suffix
    drc_match = dir3+'drc_match_'+suffix

    filters = ['F606W','F814W']

    my_header = 'RA_f606w DEC_f606w flux_f606w c_star_f606w x_f606w y_f606w rawMean_f606w zptMean_f606w rawErr_f606w zptErr_f606w  pos_std_f606w RA_f814w DEC_f814w flux_f814w c_star_f814w x_f814w y_f814w rawMean_f814w zptMean_f814w rawErr_f814w zptErr_f814w  pos_std_f814w'

    # my_form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f'

    drc_header = "RA_v DEC_v x_v y_v fAper_v fErr_v magAper_v magErr_v magRaw_v magRed_v magAbs_v elong_v ellip_v class_Star_v RA_i DEC_i x_i y_i fAper_i fErr_i magAper_i magErr_i magRaw_i magRed_i magAbs_i elong_i ellip_i class_Star_i corrF_errV corrF_errI corrM_errV corrM_errI"

    for ff in range(len(filters)):
        dataFile = np.genfromtxt(my_file)
        idxFile = np.genfromtxt(my_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir3+'my'+suffix,outArr,header=my_header)

        dataFile = np.genfromtxt(drc_file)
        idxFile = np.genfromtxt(drc_match,dtype=int)

        outArr = dataFile[idxFile]

        np.savetxt(dir3+'drc'+suffix,outArr,header=drc_header)

    return None

run_sph(drc_file,my_file,suffix='_drcLast.dat')
matchMags('_drcLast.dat')
# run_sph(drc_file,my_file,ra_idx=-2,dec_idx=-1,suffix='Dith234.dat')
# matchMags('Dith234.dat')
