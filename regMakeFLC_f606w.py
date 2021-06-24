import numpy as np

workDir='./catRawMags1305/catDir/'

file = np.genfromtxt(workDir+'matched_w_MagsPos1106r2.dat')

RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4, id = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41

newCol=np.zeros((len(file),1))
newCol[:,0] = np.arange(0,len(file),1)

file = np.hstack((file,newCol))

idx= np.argsort(file[:,mag1])[:100]

file_100 = file[idx]

np.savetxt('flc100source_f606w.reg',file_100[:,[RA,DEC]],fmt='%1.7f')

# header = 'RA DEC flux flags c_star mag1 mag2 mag3 mag4 ra1 '
# header += 'dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xr2 yr2 xr3 yr3 xr4 yr4 '
# header += 'xc1 yc1 xc2 yc2 xc3 yc3 xc4 yc4 xt1 yt1 xt2 yt2 xt3 yt3 xt4 yt4'
#
# form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f '
# form +='%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f '
# form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
# form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
# form +='%1.4f %1.4f'

# np.savetxt('flc6source.dat',file_50[[8,9,10,12,16,18],:],fmt=form,header=header)

np.savetxt(workDir+'flc100_f606w.dat',file_100[:,[RA,DEC,xt1,yt1,id]],fmt='%1.7f  %1.7f  %1.5f %1.5f %d')
