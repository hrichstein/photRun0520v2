import numpy as np

matchtol = 0.05

dir = 'sDRC_2606/'

f606w = np.genfromtxt(dir+'cat_HOR-I_F606W_cStarCut_bkgd.dat')
f814w = np.genfromtxt(dir+'cat_HOR-I_F814W_cStarCut_bkgd.dat')


flags, RA, DEC, xr, yr, flux, bkgd, c_star, magr, id, idNew = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10


f606w_id = np.zeros((len(f606w),1))
f606w_id[:,0] = np.arange(0,len(f606w),1) # index, DRCx DRCy magDRC
f606w = np.hstack((f606w,f606w_id))

f814w_id = np.zeros((len(f814w),1))
f814w_id[:,0] = np.arange(0,len(f814w),1) # index, DRCx DRCy magDRC
f814w = np.hstack((f814w,f814w_id))

master = f814w[:,[xr,yr,magr,idNew]]

x,y,magr = 0,1,2

matchids = np.zeros((len(master),1)) #id_f606w
#xDRC yDRC magDRC


nF = True
row = 0

while (nF):
    matchrows = f606w[(abs(master[row][x] - f606w[:,xr]) <= matchtol) & (abs(master[row][y] - f606w[:,yr]) <= matchtol)]

    if (len(matchrows) == 1):
      matchids[row][0] = matchrows[0][idNew]
      row = row + 1

    # elif (len(matchrows) > 1):
    #      distDiff = np.zeros((len(matchrows),1))
    #      for dd in range(len(matchrows)):
    #          distDiff[dd] = np.sqrt( (master[row][x] - matchrows[dd][xr])**2 +  (master[row][y] - matchrows[dd][yr])**2)
    #      small = np.argmin(distDiff)
    #      matchids[row][0] = matchrows[small][id]
    #      row += 1

    else:
        master = np.delete(master,row,0)
        matchids = np.delete(matchids,row,0)

    if (row >= len(master)):
        nF = False
        print(len(master))


master= np.hstack((master,matchids))
header= 'x_f814w y_f814w magr_f814w id_f814w id_f606w'


np.savetxt(dir+'drcRefidx_2906_all_cut_bkgd.dat',master,header=header)
