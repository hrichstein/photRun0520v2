import numpy as np
from linear6d import *
from getJdan import *

dir = 'catRawMags1305/catDir/'
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


dir3  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
drc_file = np.genfromtxt(dir3 + 'HOROLOGIUM-I_sfErr.dat')

RA_v, DEC_v, x_v, y_v, fAper_v, fErr_v, magAper_v, magErr_v, magRaw_v, magRed_v, magAbs_v, elong_v, ellip_v, class_Star_v, RA_i, DEC_i, x_i, y_i, fAper_i, fErr_i, magAper_i, magErr_i, magRaw_i, magRed_i, magAbs_i, elong_i, ellip_i, class_Star_i, corrF_errV, corrF_errI, corrM_errV, corrM_errI,id = 0, 1, 2 ,3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32


x = drc_file[:,x_v]
y = drc_file[:,y_v]
# Original DRC X,Y
x_arr = np.array([x[834],x[523],x[338],x[467],x[81],x[44],x[0]])
y_arr = np.array([y[834],y[523],y[338],y[467],y[81],y[44],y[0]])


match_arr = np.zeros((len(x_arr),2))

match_arr[:,0] = x_arr
match_arr[:,1] = y_arr


# PSF positions (subbing for aper)
x_mas = np.array([4718.1953,3991.6089,3783.894,3711.5825,2812.5161,2716.2554,2434.9197 ])

y_mas = np.array([4569.6943,5523.9956,5660.8569,5374.0845,5615.8149,5347.2842,5232.4946])

master_arr = np.zeros((len(x_mas),2))

master_arr[:,0] = x_mas
master_arr[:,1] = y_mas

def linTrans(targname,filt,dir='catRawMags1305/catDir/'):

    match = match_arr
    master= master_arr

    weights = np.zeros((len(match)))
    weights.fill(1.0)


    all = drc_file


    outname = "drc2APER_f606w.dat"

    new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,x_v],all[:,y_v])

    np.savetxt(dir+outname, new_all, fmt="%1.6f")


    return None

def addTranscols(targname,filt,dir='catRawMags1305/catDir/'):

    cat = drc_file
#
    transCat = np.loadtxt(dir+"drc2APER_f606w.dat")
#
    newCol = np.zeros((len(cat),2))

    newCol[:,0] = transCat[:,0]
    newCol[:,1] = transCat[:,1]

    cat = np.hstack((cat, newCol))
#
    header = 'RA_v DEC_v x_v y_v fAper_v fErr_v magAper_v magErr_v magRaw_v magRed_v magAbs_v elong_v ellip_v class_Star_v RA_i DEC_i x_i y_i fAper_i fErr_i magAper_i magErr_i magRaw_i magRed_i magAbs_i elong_i ellip_i class_Star_i corrF_errV corrF_errI corrM_errV corrM_errI xt yt'

    np.savetxt(dir+'drc2APER_f606w_t.dat',cat,header=header)
#
#
    return None

targname = 'HOROLOGIUM-I'
filt = 'F606W'

linTrans(targname,filt,dir='catRawMags1305/catDir/')
addTranscols(targname,filt,dir='catRawMags1305/catDir/')
