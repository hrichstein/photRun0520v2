import numpy as np

flags_f606w, RA_f606w, DEC_f606w, xr_f606w, yr_f606w, flux_f606w, bkgd_f606w, c_star_f606w, magr_f606w, id_f606w, xr_f814w_trans, yr_f814w_trans, flags_f814w, RA_f814w, DEC_f814w, xr_f814w, yr_f814w, flux_f814w, bkgd_f814w, c_star_f814w, magr_f814w, id_f814w,id = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22

dir = 'sDRC_2606/'

drcCat = np.genfromtxt(dir +'matchedDRCfullCat_2906_bkgd.dat')

sort_idx = np.argsort(drcCat[:,magr_f606w])[:50]

idxCol = np.zeros((len(sort_idx),1))
idxCol[:,0] = np.arange(0,len(sort_idx),1)

drc = drcCat[sort_idx]
drc= np.hstack((drc,idxCol))

np.savetxt('drc50_bkgd.reg',drc[:,[RA_f606w,DEC_f606w]],fmt='%1.7f')

np.savetxt('drc50_bkgd.dat',drc[:,[RA_f606w,DEC_f606w,xr_f606w,yr_f606w,id]],fmt='%1.7f %1.7f %1.5f %1.5f %d')
