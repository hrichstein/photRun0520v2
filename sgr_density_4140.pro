PRO sgr_density
;
; Distance to the Sgr stars
;
dist = 33.4*1E3
mu = 5.0D * ALOG10(dist) - 5.0D ;dist modulus
;
; Read in observed data and isochrone.
; "sgr1.dat" was queried using the following SQL command on SkyServer DR10:
; (http://skyserver.sdss3.org/dr10/en/tools/search/sql.aspx)
;
; SELECT s.ra, s.dec, s.psfmag_g-s.extinction_g as g, 
;        s.psfmag_r-s.extinction_r as r, s.psfmag_i-s.extinction_i as i
;        FROM star s
;        WHERE 
;        s.ra BETWEEN  31.9 AND 33.0 and s.dec BETWEEN -5.0 and -4.0
;        AND CLEAN=1 AND s.psfmag_r <23
;
readcol,'sgr1.dat',f='d,d,d,d,d',ra,dec,gobs,robs,iobs,/silent ;;; De-reddened magnitudes
readcol,'sdss.iso',f='i,d,d,d,d,d,d,d,d,d',no,no,no,no,no,uiso,giso,riso,iiso,ziso,/silent
giso=giso+mu
riso=riso+mu
iiso=iiso+mu
;
; Select stars in SDSS catalog that are consistent with isochrone population at 33.4 kpc.
; Get everything within 0.1 mags of above isochrone and also where sdss has 95% completeness
;
distmin=DBLARR(n_elements(iobs))
for i = 0,n_elements(iobs)-1 do begin
   dist = SQRT((iobs[i]-iiso)^2 + ((gobs[i]-iobs[i])-(giso-iiso))^2)
   distmin[i]= min(dist)
endfor
s=where(distmin LE 0.1D and iobs LE 21.3D,sdss_n)
;
; According to website below, DR10 data is 95% complete down to g=22.2, r=22.2, i=21.3
; http://www.sdss3.org/dr10/scope.php#opticalstats
; Therefore, I chose i mag down to 21.3.
;;;;
;
; Now the background. For the background, I chose the same Galactic longitude, but opposite latitude.
;
readcol,'background.dat',f='d,d,d,d,d',bra,bdec,bgobs,brobs,biobs,/silent
bdistmin=DBLARR(n_elements(biobs))
for i = 0,n_elements(biobs)-1 do begin
   bdist = SQRT((biobs[i]-iiso)^2 + ((bgobs[i]-biobs[i])-(giso-iiso))^2)
   bdistmin[i]= min(bdist)
endfor
bs=where(bdistmin LE 0.1D and biobs LE 21.3D,bsdss_n)
;;print,sdss_n
;;print,bsdss_n
;
; Do the background subtraction
;
sdss_n = sdss_n-bsdss_n
;
; Now "sdss_n" is the predictied number of Sgr stream like population after background subtraction.
;
;;;;
;
; We now have to use theoretical LFs to scale the number of stars in i<21.3
; down to whichever magnitude range consistent with our HST mags (20.0<F814W<25.5).
;
; From Sirianni et al. (2005) with a little basic math for FIELD1:
; I = (F814W-25.531)+25.493 -0.003(F606W-F814W)
; V = (F606W-26.407)+26.582 +0.310(F606W-F814W)
;
; Webpage http://www.astro.ex.ac.uk/people/timn/Catalogues/Sloan_vs_Cousins/ gives:
; V-i = 0.971(V-I) + 0.008
; so, i = V -0.971(V-I) -0.008
; This might not be that accurate, but it'll be enough for our purpose.
; Using equations above give our pm_field1_groupab.txt i_sdss range of:
; 20.9 < i_sdss < 25.5
;
;;;
;;; Deal with LF
;;;
readcol,'SDSS.lf.i.txt',f='i,d,d,d',no,iband,logn,logdn,/silent
N=10^(logn)
;
; N is cumulative number (that is, number of stars brighter than iband[i])

iband=iband+mu
s=where(iband lt 21.1 and iband gt 20.9)
ns = N[s[0]]

;Now calculate number of stars in desired region
s1=where(iband lt 25.6 and iband gt 25.4)
s2=where(iband lt 20.9 and iband gt 20.7)
;;print,N[s1[0]]
;;print,N[s2[1]]
ns2=N[s1[0]]-N[s2[1]] ; Number in our observed mag region

; ns:ns2 = sdss_n:X
; X = (ns2*sdss_n)/ns

WFC_FOV=(202.^2)/(3600.^2) ;degrees
print,(ns2*sdss_n)/ns*WFC_FOV)
;
END
