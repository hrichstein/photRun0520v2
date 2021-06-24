import numpy as np

x = np.arange(0,1000,1)
y = np.arange(0,2000,1)

good = (x >= 0) & (y >= 0)

# Best fit parameters
richness = 400
x0 = 4000
y0 = 3600
gal_ext_pix = 5200
gal_ell = 0.5
gal_pa = 80  # in degrees

# Define pixel grid

# 267 is number of bins, I think
xbin0 = np.ones(267)
xbin1 = np.arange(0,267,1)
xbin2 = xbin1 * 30 + 15

xbin = np.matmul(xbin0,xbin2)
ybin = np.matmul(xbin2,xbin0)

# plummer model
costh = np.cos(-1 * np.deg2rad(gal_pa))
sinth = np.sin(-1 * np.deg2rad(gal_pa))
dx = xbin - x0
dy = ybin - y0

r1 = (dx * costh - dy * sinth) / (1 - gal_ell)
r2 = (dx * sinth + dy * costh)

radius = np.sqrt(r1**2 + r2**2)

r_h = gal_ext_pix
norm = r_h**2 / (np.pi * (1-gal_ell))
pdf = norm/((radius**2 + r_h**2)**2)

xdel = 30   # Figure out where this 30 comes from....
ydel = 30
pixarea = xdel * ydel

model_counts_gal = richness * pdf * pixarea

# Observed surface density map
binned_surface_density_30 = np.zeros((267,267))
for ii in range(267):
    for jj in range(267):
        idx = (x[good] >= 30 * ii & x[good] <= 30 * (ii + 1)) \
            & (y[good] >= 30 * jj & y[good] <= 30 * (jj + 1))

        tmp = np.where(idx)[0]
        # flipping ii and jj because in IDL, it's col, row.
        binned_surface_density_30[jj,ii] = len(tmp)

# Define annuli
r1 = (radius < (0.1 * r_h)).nonzero()
r2 = ((radius < (0.2 * r_h)) & (radius >= (0.1 * r_h))).nonzero()
r3 = ((radius < (0.3 * r_h)) & (radius >= (0.2 * r_h))).nonzero()
r4 = ((radius < (0.4 * r_h)) & (radius >= (0.3 * r_h))).nonzero()
r5 = ((radius < (0.5 * r_h)) & (radius >= (0.4 * r_h))).nonzero()
r6 = ((radius < (0.6 * r_h)) & (radius >= (0.5 * r_h))).nonzero()
r7 = ((radius < (0.7 * r_h)) & (radius >= (0.6 * r_h))).nonzero()
r8 = ((radius < (0.8 * r_h)) & (radius >= (0.7 * r_h))).nonzero()
r9 = ((radius < (0.9 * r_h)) & (radius >= (0.8 * r_h))).nonzero()
r10 = ((radius < r_h) & (radius >= (0.9 * r_h))).nonzero()

ra1 = np.arange(0,10,1)
r_annuli = (ra1 + 0.05) * r_h

area_annuli = np.pi * ((r_annuli + 0.05)**2 - (r_annuli - 0.05)**2)

data_1d = [binned_surface_density_30[r1].sum(),
           binned_surface_density_30[r2].sum(),
           binned_surface_density_30[r3].sum(),
           binned_surface_density_30[r4].sum(),
           binned_surface_density_30[r5].sum(),
           binned_surface_density_30[r6].sum(),
           binned_surface_density_30[r7].sum(),
           binned_surface_density_30[r8].sum(),
           binned_surface_density_30[r9].sum(),
           binned_surface_density_30[r10].sum()]

model_1d = [model_counts_gal[r1].sum(),
            model_counts_gal[r2].sum(),
            model_counts_gal[r3].sum(),
            model_counts_gal[r4].sum(),
            model_counts_gal[r5].sum(),
            model_counts_gal[r6].sum(),
            model_counts_gal[r7].sum(),
            model_counts_gal[r8].sum(),
            model_counts_gal[r9].sum(),
            model_counts_gal[r10].sum()]






















#
