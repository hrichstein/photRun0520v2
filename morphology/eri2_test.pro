fname = '/Volumes/Spare Data/photRun0520/morphology/eri2cut.cat'

; OPENR,1,fname
ASTR=READ_ASCII(fname,DATA_START=1,HEADER=h)
A=ASTR.(0)
x=A[0,*]
y=A[1,*]
m606=A[2,*]
m814=A[4,*]
flag=A[6,*]

;BEST-FIT ERI II PARAMETERS (FROM MCMC_STRUCTURAL_FIT_ERI2_FINAL.IPYNB)
richness = 29662
x0 = 4123.4
y0 = 3601.7
eri_ext_pix = 5205.
eri_ell = 0.447
eri_pa = 77.8

;DFINE PIXEL GRID
xbin = (fltarr(267)+1)##(findgen(267)*30. + 15)
; rows go from 15 to 315+, left to right
; columns down are the same as the number in the first row
; a b c
; a b c
; a b c

ybin = (findgen(267)*30.+15)##(fltarr(267)+1)
; rows are the same from left to right
; columns increase by 30
; a a a
; b b b
; c c c

;PLUMMER MODEL OF ERI II
costh = cos(-1*eri_pa*!dtor)  ;dtor is degrees to radians conversion
sinth = sin(-1*eri_pa*!dtor)
dx = xbin-x0  ;the differences between each xbin center and the x-origin
dy = ybin-y0  ;the differences between each ybin center and the y-origin
radius = sqrt( ((dx*costh - dy*sinth)/(1-eri_ell))^2 + (dx*sinth + dy*costh)^2)

r_h = eri_ext_pix
norm = r_h^2/(!pi*(1 - eri_ell))
pdf = norm/((radius^2 + r_h^2)^2)

xdel = 30.
ydel = 30.
pixarea = xdel*ydel

model_counts_eri2 = richness*pdf*pixarea


;OBSERVED SURFACE DENSITY MAP
binned_surface_density_30 = fltarr(267,267)
for i = 0,266 do begin &$
   for j = 0,266 do begin &$
     inbin = where(x ge 30*i AND x lt 30*(i+1) AND $
                   y ge 30*j AND y lt 30*(j+1),n_inbin) &$
     binned_surface_density_30[i,j] = n_inbin &$
   endfor &$
endfor


;DEFINE ANNULI - The where looks at radius like a flattened array
r1 = where(radius lt 0.1*r_h)
r2 = where(radius lt 0.2*r_h AND radius ge 0.1*r_h)
r3 = where(radius lt 0.3*r_h AND radius ge 0.2*r_h)
r4 = where(radius lt 0.4*r_h AND radius ge 0.3*r_h)
r5 = where(radius lt 0.5*r_h AND radius ge 0.4*r_h)
r6 = where(radius lt 0.6*r_h AND radius ge 0.5*r_h)
r7 = where(radius lt 0.7*r_h AND radius ge 0.6*r_h)
r8 = where(radius lt 0.8*r_h AND radius ge 0.7*r_h)
r9 = where(radius lt 0.9*r_h AND radius ge 0.8*r_h)
r10 = where(radius lt 1*r_h AND radius ge 0.9*r_h)

r_annuli = (findgen(10) + 0.05)*r_h
area_annuli = !pi*( (r_annuli + 0.05)^2 - (r_annuli - 0.05)^2)
data_1d = [total(binned_surface_density_30[r1]), $
            total(binned_surface_density_30[r2]), $
            total(binned_surface_density_30[r3]), $
            total(binned_surface_density_30[r4]), $
            total(binned_surface_density_30[r5]), $
            total(binned_surface_density_30[r6]), $
            total(binned_surface_density_30[r7]), $
            total(binned_surface_density_30[r8]), $
            total(binned_surface_density_30[r9]), $
            total(binned_surface_density_30[r10])]
model_1d = [total(model_counts_eri2[r1]), $
             total(model_counts_eri2[r2]), $
             total(model_counts_eri2[r3]), $
             total(model_counts_eri2[r4]), $
             total(model_counts_eri2[r5]), $
             total(model_counts_eri2[r6]), $
             total(model_counts_eri2[r7]), $
             total(model_counts_eri2[r8]), $
             total(model_counts_eri2[r9]), $
             total(model_counts_eri2[r10])]

end
