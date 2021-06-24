import numpy as np

dir = 'catRawMags1305/catDir/'

# master= np.genfromtxt(dir+'psfPSF_idx_1106_mas.dat')
master= np.genfromtxt(dir+'drcAPER_idx_2406_mas.dat')


f606w = np.genfromtxt(dir+'matchedDRCaper2406_606.dat') # 606
f814w = np.genfromtxt(dir+'matchedDRCaper2406_814.dat') #814

# There should be a last column in each catalog giving the index placement

mag606, magErr606, xt_v, yt_v, xAPER_606, yAPER_606, magAPER_606, idAPER_606,id_new606 = 0, 1, 2, 3, 4, 5, 6, 7,8

mag814, magErr814, xt_v, yt_v, xAPER_814, yAPER_814, magAPER_814, idAPER_814,id_new814 = 0, 1, 2, 3, 4, 5, 6, 7,8

x,y,mag606_in,id606,id814 = 0,1,2,3,4
# Master header
# x y magZ id606 id814

# newCols606 = np.zeros((len(master),11))

idCol606 = master[:,id606]
idx606 = np.asarray(idCol606,int)
reg606 = f606w[idx606]

##################################

# newCols814 = np.zeros((len(master),11))
idCol814 = master[:,id814]
idx814 = np.asarray(idCol814,int)
reg814 = f814w[idx814]

outArr = np.hstack((reg606,reg814))
header= 'mag606 magErr606 xt_v yt_v xAPER_606 yAPER_606 magAPER_606, idAPER_606 mag814 magErr814 xt_i yt_i xAPER_814 yAPER_814 magAPER_814, idAPER_814 '


# np.savetxt(dir+'matchedPSFpsf1106.dat',outArr,header=header)
np.savetxt(dir+'matchedAPERdrc2406_comb.dat',outArr,header=header)
