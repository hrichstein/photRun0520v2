import numpy as np

dir = 'catRawMags1305/catDir/'

# 0506 run
# data = np.genfromtxt(dir+'flcDRCpos0506.dat')
# 1106 run
data = np.genfromtxt(dir+'flcDRCpos1106_f606wr2.dat')

stdCut=2.5
RA, DEC, flux, flags, c_star, mag1, mag2, mag3, mag4, ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4, xr1, yr1, xr2, yr2, xr3, yr3, xr4, yr4, xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, xt1, yt1, xt2, yt2, xt3, yt3, xt4, yt4, xDRC, yDRC = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42

newCol = np.zeros((len(data),5))

counter = 0
keptAll = 0
naSTD = int(0)
n3STD = int(0)
for dd in range(len(data)):
    idx = np.array([mag1,mag2,mag3,mag4],dtype=int) # columns of mags

    tmp_std_arr = np.zeros((4))

    arr1 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[2]]]) #3 is out
    arr2 = np.array([data[dd,idx[3]], data[dd,idx[1]], data[dd,idx[2]]]) #0 is out
    arr3 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[3]]]) #2 is out
    arr4 = np.array([data[dd,idx[0]], data[dd,idx[2]], data[dd,idx[3]]]) #1 is out

    mean1 = np.nanmean(arr1)
    mean2 = np.nanmean(arr2)
    mean3 = np.nanmean(arr3)
    mean4 = np.nanmean(arr4)

    std1 = np.nanstd(arr1)
    std2 = np.nanstd(arr2)
    std3 = np.nanstd(arr3)
    std4 = np.nanstd(arr4)

    tmp_std_arr[0] = abs(data[dd,idx[3]] - mean1)/std1
    tmp_std_arr[1] = abs(data[dd,idx[0]] - mean2)/std2
    tmp_std_arr[2] = abs(data[dd,idx[2]] - mean3)/std3
    tmp_std_arr[3] = abs(data[dd,idx[1]] - mean4)/std4

    max_std = np.nanmax(tmp_std_arr)

    num_above_std = (tmp_std_arr >= stdCut).sum()
    num_above_std = int(num_above_std)

    if num_above_std == 2:
        # print('Num Above Std: ', num_above_std)
        naSTD += 1
    if num_above_std > 2:
        n3STD += 1

    if max_std >= stdCut:
        tmp_idx = np.argsort(tmp_std_arr)[::-1][0]
        if tmp_idx == 0:
            newCol[dd,0] = mean1
            newCol[dd,1] = std1
            newCol[dd,2] = 1 # flag cut
            newCol[dd,3] = 3 # idx cut
            counter += 1

        elif tmp_idx == 1:
            newCol[dd,0] = mean2
            newCol[dd,1] = std2
            newCol[dd,2] = 1
            newCol[dd,3] = 0
            counter += 1

        elif tmp_idx == 2:
            newCol[dd,0] = mean3
            newCol[dd,1] = std3
            newCol[dd,2] = 1
            newCol[dd,3] = 2
            counter += 1

        else:
            newCol[dd,0] = mean4
            newCol[dd,1] = std4
            newCol[dd,2] = 1
            newCol[dd,3] = 1
            counter += 1

    else:
        newCol[dd,0] = np.nanmean(data[dd][idx])
        newCol[dd,1] = np.nanstd(data[dd][idx])
        newCol[dd,2] = 0
        newCol[dd,3] = 4 # means no index was cut
        keptAll += 1

    newCol[dd,4] = num_above_std

out = data[:,[RA, DEC, flags, c_star, mag1, mag2, mag3, mag4, xr1, yr1, xt1, yt1]]

data = np.hstack((out, newCol))

header = 'RA DEC flags c_star mag1 mag2 mag3 mag4 xr1 yr1 xt1 yt1 mean stdev cut_flag idx_cut num_abv_std'

# np.savetxt('flcAll_magCut_0906.dat',data,header=header)
np.savetxt('flcAll_magCut_1106.dat',data,header=header)

print('Counter:',counter)
print('Same: ',keptAll)
print('Num with 2 above STD: ', naSTD)
