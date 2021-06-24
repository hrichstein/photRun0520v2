import numpy as np
from linear6d import *
from getJdan import *

dir = 'sDRC_2606/'

dir2 = 'catRawMags1305/catDir/'
'''

transform.py

This code serves as the main interface for using the mpfit Python code (mpfit.py) to
determine the coefficients for a 6-parameter linear transformation (linear6d.py).
For the input, it takes 8 different 1D arrays:
- x and y positions of a list of sources to be transformed into another frame that
have already been matched to a corresponding list of sources in the other frame.

- x and y positions of the corresponding list in the intended frame

- weights for the sources (this is currently set to 1.0 for every source
to give equal weight

- x and y positions of all sources that are to be transformed into another frame. Once
the 6-param coefficients have been determined using the above matched sources, the
transformation is applied to this full list of sources and is saved out.


The first set of coordinates is the reference sources in your particular secondary dither. The second set of coordinates is the reference sources in the primary dither. And the third set is all sources in your secondary dither.
'''

# APER

# xPSF, yPSF, m606cPSF, m814cPSF, s606PSF, s814PSF, nstarPSF, nstarAPER, idAPER, xAPER, yAPER, m606cAPER, m814cAPER, s606APER, s814APER =  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14

# aperCat = np.genfromtxt(dir2+'matchedPSFaper2906_tc.dat')
aperReg = np.genfromtxt('aper50_bkgd.dat')

aperArr = np.array([aperReg[2],aperReg[3],aperReg[6],aperReg[7],aperReg[10],aperReg[18],aperReg[21]])

master_arr = np.zeros((len(aperArr),2))

master_arr[:,0] = aperArr[:,2]
master_arr[:,1] = aperArr[:,3]
# DRC

flags_f606w, RA_f606w, DEC_f606w, xr_f606w, yr_f606w, flux_f606w, bkgd_f606w, c_star_f606w, magr_f606w, id_f606w, xr_f814w_trans, yr_f814w_trans, flags_f814w, RA_f814w, DEC_f814w, xr_f814w, yr_f814w, flux_f814w, bkgd_f814w, c_star_f814w, magr_f814w, id_f814w = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21

# 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40

drcCat = np.genfromtxt(dir +'matchedDRCfullCat_2906_bkgd.dat')

drcReg = np.genfromtxt('drc50_bkgd.dat')

drcArr = np.array([drcReg[23],drcReg[25],drcReg[27],drcReg[28],drcReg[32],drcReg[39],drcReg[41]])

match_arr = np.zeros((len(drcArr),2))

match_arr[:,0] = drcArr[:,2]
match_arr[:,1] = drcArr[:,3]

print(len(drcArr),len(aperArr))

def linTrans(targname,filt):

    match = match_arr
    master= master_arr

    weights = np.zeros((len(match)))
    weights.fill(1.0)

    all = drcCat

    outname = 'drc2APERmatch2906_full_bkgd.dat'

    new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,xr_f606w],all[:,yr_f606w])

    np.savetxt(dir+outname, new_all, fmt="%1.6f")
    # np.savetxt(dir+'flcDRC0406newMatch.dat',new_match,fmt='%1.6f')


    return None

def addTranscols(targname,filt):

    # cat = np.genfromtxt(dir+'matched_w_MagsPos1106r2.dat')
    # transCat = np.genfromtxt(dir + 'flcDRCmatch2506_f606w.dat')

    cat = np.genfromtxt(dir+'matchedDRCfullCat_2906_bkgd.dat')
    transCat = np.genfromtxt(dir + 'drc2APERmatch2906_full_bkgd.dat')

    newCol = np.zeros((len(cat),2))

    newCol[:,0] = transCat[:,0]
    newCol[:,1] = transCat[:,1]

    cat = np.hstack((cat, newCol))

    header = 'flags_f606w RA_f606w DEC_f606w xr_f606w yr_f606w flux_f606w bkgd_f606w c_star_f606w magr_f606w id_f606w xr_f814w_trans yr_f814w_trans flags_f814w RA_f814w DEC_f814w xr_f814w yr_f814w flux_f814w bkgd_f814w c_star_f814w magr_f814w id_f814w xAPER_trans yAPER_trans'

    # np.savetxt(dir+'flcDRCpos2506_f606w.dat', cat, header=header)

    np.savetxt(dir+'aperDRCpos2906_fullTransAPER_bkgd.dat', cat, header=header)
    return None

targname = 'HOROLOGIUM-I'
filt = 'F814W'

linTrans(targname,filt)
addTranscols(targname,filt)

# Former
# Iter      23    CHI-SQUARE =  21.01153291  DOF =  12
#    P0 = 1842.90129
#    P1 = 0.09289828656
#    P2 = 0.9950916764
#    P3 = 7105.368901
#    P4 = -0.9954595125
#    P5 = 0.09355332827

#New
# Iter      22    CHI-SQUARE =  12.22172326  DOF =  8
#    P0 = 1843.814373
#    P1 = 0.09250721194
#    P2 = 0.9948102152
#    P3 = 7104.728628
#    P4 = -0.995265243
#    P5 = 0.09381963651
