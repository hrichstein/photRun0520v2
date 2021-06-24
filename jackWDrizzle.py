from stwcs import updatewcs
import glob
import os

_ = os.system ("crds bestrefs --update-bestrefs --sync-references=1 --files *fl?.fits")
input_images = sorted(glob.glob('*fl?.fits'))
_ = updatewcs.updatewcs (input_images)
del updatewcs

from astropy.io import fits
import numpy as np

exptime = np.zeros(len(input_images))
for i in range(len(input_images)):
    fname = input_images[i]
    hdu = fits.open (fname)
    exptime[i] = hdu[0].header['expend'] - hdu[0].header['expstart']
    hdu.close()
del fits

fname_minexp = input_images[exptime.argmin()]
fname_maxexp = input_images[exptime.argmax()]

_ = os.system ("python3 hgwell.py %s 40000 30 0.5" %fname_minexp)
_ = os.system ("python3 hgwell.py %s 40000 30 0.5" %fname_maxexp)

dfs = [pd.read_csv(j) for j in sorted(glob.glob("*SCI?.csv"))]
for i,df in enumerate(dfs):
    df = df[df.g.notna()].reset_index(drop=True)
    df['sep'] = np.sqrt( (df.ra_pix - df.ra_prop)**2 + (df.dec_pix - df.dec_prop)**2 )
    dfs[i] = df
dfs = pd.concat (dfs, ignore_index=True)
dfs = dfs.sort_values (by='sep', ascending=True)
dfs = dfs.drop_duplicates (subset='g', keep='first', ignore_index=True)
dfs.to_csv ('gaia_full_match.csv', index=False)
from astropy.table import Table
tbl = Table([dfs.ra_prop, dfs.dec_prop])
tbl.write ('gaia.cat', format='ascii.fast_commented_header')
del Table,dfs,tbl

from drizzlepac import tweakreg
from stsci.tools import teal

cat = 'gaia.cat'
wcsname = 'GAIA'
teal.unlearn('tweakreg')
teal.unlearn('imagefindpars')
cw = 3.5
tweakreg.TweakReg(input_images, # Pass input images
                  updatehdr=True, # update header with new WCS solution
                  imagefindcfg={'threshold':250.,'conv_width':cw},# Detection parameters, threshold varies for different data
                  separation=0.0, # Allow for very small shifts
                  refcat=cat, # Use user supplied catalog (Gaia)
                  clean=True, # Get rid of intermediate files
                  interactive=False,
                  see2dplot=False,
                  shiftfile=True, # Save out shift file (so we can look at shifts later)
                  wcsname=wcsname, # Give our WCS a new name
                  reusename=True,
                  fitgeometry='general') # Use the 6 parameter fit
