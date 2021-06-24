import numpy as np

dir = 'sDRC_2606/'

f606w = np.genfromtxt(dir+'cat_HOROLOGIUM-I_F606W_wMag.dat')
f814w = np.genfromtxt(dir+'cat_HOROLOGIUM-I_F814W_wMag.dat')

flags, RA, DEC, xr, yr, flux, c_star, magr, id = 0, 1, 2, 3, 4, 5, 6, 7, 8

f606w_cStar = f606w[:,c_star]
f814w_cStar = f814w[:,c_star]

carry_f606w = f606w[f606w_cStar>=0.5]
carry_f814w = f814w[f814w_cStar>=0.5]

print(len(carry_f606w))
print(len(carry_f814w))

header = "flags RA DEC xr yr flux c_star magr id"

np.savetxt(dir+'cat_HOR-I_F606W_cStarCut.dat',carry_f606w,fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d",header=header)

np.savetxt(dir+'cat_HOR-I_F814W_cStarCut.dat',carry_f814w,fmt="%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d",header=header)
