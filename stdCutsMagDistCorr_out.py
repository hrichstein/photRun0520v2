import numpy as np
from hst_fun_2 import *
import matplotlib.pyplot as plt


def magCorr(dataN,filt):

    fileDir = '../Hannah_Data/'
    cor = hst_gc_z(fileDir+'wfc_'+filt,dataN['xr1'],dataN['yr1'])


    return cor


def makeSTDcuts(dir,targname,filt,suffix='_aftLT.dat'):

    dataN = np.genfromtxt(dir+'matched_w_MagsPos_'+filt+suffix,names=True)
    data = np.genfromtxt(dir+'matched_w_MagsPos_'+filt+suffix)

    colNs = np.array(dataN.dtype.names)

    mag1 = np.int(np.where(colNs=='mag1')[0])
    mag2 = np.int(np.where(colNs=='mag2')[0])
    mag3 = np.int(np.where(colNs=='mag3')[0])
    mag4 = np.int(np.where(colNs=='mag4')[0])

    stdCut=3
    cor = magCorr(dataN,filt)

    # newCol = np.zeros((len(data),5))
    #
    # counter = 0
    # keptAll = 0
    # naSTD = int(0)
    #
    # for dd in range(len(data)):
    #     idx = np.array([mag1,mag2,mag3,mag4],dtype=int) # columns of mags
    #     tmp_std_arr = np.zeros((4))
    #
    #     arr1 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[2]]]) #3 is out
    #     arr2 = np.array([data[dd,idx[3]], data[dd,idx[1]], data[dd,idx[2]]]) #0 is out
    #     arr3 = np.array([data[dd,idx[0]], data[dd,idx[1]], data[dd,idx[3]]]) #2 is out
    #     arr4 = np.array([data[dd,idx[0]], data[dd,idx[2]], data[dd,idx[3]]]) #1 is out
    #
    #     mean1 = np.nanmean(arr1) + cor[0]
    #     mean2 = np.nanmean(arr2) + cor[0]
    #     mean3 = np.nanmean(arr3) + cor[0]
    #     mean4 = np.nanmean(arr4) + cor[0]
    #
    #     std1 = np.nanstd(arr1)
    #     std2 = np.nanstd(arr2)
    #     std3 = np.nanstd(arr3)
    #     std4 = np.nanstd(arr4)
    #
    #     tmp_std_arr[0] = abs(data[dd,idx[3]] - mean1)/std1
    #     tmp_std_arr[1] = abs(data[dd,idx[0]] - mean2)/std2
    #     tmp_std_arr[2] = abs(data[dd,idx[2]] - mean3)/std3
    #     tmp_std_arr[3] = abs(data[dd,idx[1]] - mean4)/std4
    #
    #     max_std = np.nanmax(tmp_std_arr)
    #
    #     num_above_std = (tmp_std_arr >= stdCut).sum()
    #     num_above_std = int(num_above_std)
    #
    #     if num_above_std == 2:
    #         naSTD += 1
    #
    #     if max_std >= stdCut:
    #         tmp_idx = np.argsort(tmp_std_arr)[::-1][0]
    #         if tmp_idx == 0:
    #             newCol[dd,0] = mean1
    #             newCol[dd,1] = std1
    #             newCol[dd,2] = 1 # flag cut
    #             newCol[dd,3] = 3 # idx cut
    #             counter += 1
    #
    #         elif tmp_idx == 1:
    #             newCol[dd,0] = mean2
    #             newCol[dd,1] = std2
    #             newCol[dd,2] = 1
    #             newCol[dd,3] = 0
    #             counter += 1
    #
    #         elif tmp_idx == 2:
    #             newCol[dd,0] = mean3
    #             newCol[dd,1] = std3
    #             newCol[dd,2] = 1
    #             newCol[dd,3] = 2
    #             counter += 1
    #
    #         else:
    #             newCol[dd,0] = mean4
    #             newCol[dd,1] = std4
    #             newCol[dd,2] = 1
    #             newCol[dd,3] = 1
    #             counter += 1
    #
    #     else:
    #         newCol[dd,0] = np.nanmean(data[dd][idx])
    #         newCol[dd,1] = np.nanstd(data[dd][idx])
    #         newCol[dd,2] = 0
    #         newCol[dd,3] = 4 # means no index was cut
    #         keptAll += 1
    #
    #     newCol[dd,4] = num_above_std
    #
    # data = np.hstack((data, newCol))
    #
    # header = 'RA DEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xr2 yr2 xr3 yr3 xr4 yr4 xc1 yc1 xc2 yc2 xc3 yc3 xc4 yc4 xt1 yt1 xt2 yt2 xt3 yt3 xt4 yt4 mean stdev cut_flag idx_cut num_abv_std'
    #
    # form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f '
    # form +='%1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f '
    # form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
    # form +='%1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f '
    # form +='%1.4f %1.4f %1.4f %1.5f %d %d %d %1.4f'
    #
    # np.savetxt(dir+'magSTDcutAll_'+filt+'_mDc_corr.dat',data,header=header,fmt=form)
    #
    # print('Counter:',counter)
    # print('Same: ',keptAll)
    # print('Num with 2 above STD: ', naSTD)
    #
    # outCor = np.zeros((len(cor),1))
    # outCor[:,0] = cor

    np.savetxt(dir+'magDistCorrs.dat',cor,fmt='%1.4e')

    return None
