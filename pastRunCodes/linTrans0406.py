import numpy as np
from linear6d import *

dir = 'catRawMags1305/catDir/'

# ref_file = np.genfromtxt(dir+'refStars0406.dat',names=True)
# ref_file = np.genfromtxt(dir+'refStars0406r2.dat',names=True)
ref_file = np.genfromtxt(dir+'refStars0406r3.dat',names=True)

# all_file = np.genfromtxt(dir+'flc2PSF_round1.dat')
# all_file = np.genfromtxt(dir+'flc2PSF_round2.dat')
all_file = np.genfromtxt(dir+'flc2PSF_round3.dat')

match_arr = np.zeros((len(ref_file),2))
master_arr = np.zeros((len(ref_file),2))

match_arr[:,0] = ref_file['x_cat']
match_arr[:,1] = ref_file['y_cat']

master_arr[:,0] = ref_file['x_psf']
master_arr[:,1] = ref_file['y_psf']

match = match_arr
master = master_arr
all = all_file

weights = np.zeros((len(match)))
weights.fill(1.0)

new_match, new_all = test_linear(match[:,0],match[:,1], master[:,0], master[:,1], weights, weights, all[:,0],all[:,1])

# outname = 'flc2PSF_round3.dat'
outname = 'flc2PSF_round4.dat'

np.savetxt(dir+outname, new_all, fmt="%1.6f")

# After round 2
# Iter      24    CHI-SQUARE =  0.7860357406  DOF =  24
#    P0 = 0.01610738823
#    P1 = 0.9999972325
#    P2 = -7.699677485e-07
#    P3 = 2.151030895
#    P4 = -0.0001302926941
#    P5 = 0.9996981841

# After round 3
# Iter      17    CHI-SQUARE =  0.7860443496  DOF =  24
#    P0 = 3.39650612e-05
#    P1 = 0.9999999931
#    P2 = -2.925792871e-10
#    P3 = -8.666368651e-05
#    P4 = 1.501954585e-08
#    P5 = 1.000000007

# After round 4
# Iter      17    CHI-SQUARE =  0.7860479692  DOF =  24
#    P0 = -4.19837242e-05
#    P1 = 1.000000014
#    P2 = -2.5754552e-09
#    P3 = -0.000277979864
#    P4 = 1.596027933e-08
#    P5 = 1.000000037
