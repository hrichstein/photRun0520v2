'''
hst_fun.py
A series of HST-focused functions to be used in other scripts.
Late Edited: 6/30/17
'''

import numpy as np
from astropy.io import fits


################################################################
#### Takes in a list/array and turns it into a single object to
#### be appended to a list. This version assumes floats.

def addnumrow(row):

###Creates the empty list
  el = []

###Cycles through the items in the array
  for m in range(len(row)):
    el.append(float(row[m]))

###Returns the new list back
  return el

################################################################
#### Converts STScI output into useable Python arrays

def stsci_convert(filename):

###Loads in the file and shifts it into an initial list
  values = np.genfromtxt(filename,dtype=None,comments="#")
  newvalues = []

###Loop through the list, checking for bad values before adding
###to new list
  for i in range(len(values)):
    if (values[i][-1] == "**********"):
      continue
    else:
      newvalues.append(addnumrow(values[i]))

### Converts the list to an array before sending it back
  newvalues = np.asarray(newvalues)

  return newvalues


################################################################

################################################################
#### Geometric distortion correction function for the WFC3UVIS Camera on HST
#### Expects a string for 'filter' containing the desired filter
#### and 'x' and 'y' should be arrays containing the float positions of the sources

def wfc3uvis_gc(filter, x, y):

###Checks for which filter is being handled and sets the appropriate files
  if (filter == 'F606W'):
    gcxname = '/Users/paulzivick/Research/HST/GC/wfc3uv_F606W_gcx.fits'
    gcyname = '/Users/paulzivick/Research/HST/GC/wfc3uv_F606W_gcy.fits'
  elif (filter == 'F814W'):
    gcxname = '/Users/paulzivick/Research/HST/GC/wfc3uv_F814W_gcx.fits'
    gcyname = '/Users/paulzivick/Research/HST/GC/wfc3uv_F814W_gcy.fits'
  
###Shifts the values to new arrays to manipulate
  xr = x
  yr = y

  ix = np.zeros(len(xr))
  iy = np.zeros(len(yr))

  for i in range(len(x)):
    ix[i] = int(xr[i])
    iy[i] = int(yr[i])

###Loads in the distortion coefficients
  gcx = fits.open(gcxname)
  gcxim = gcx[0].data
  gcx.close()

  gcy = fits.open(gcyname)
  gcyim = gcy[0].data
  gcy.close()

###Checks for oddities in the indicies of the source positions
  chip = 1
  k = ''
  
  Ns = len(xr)

  nj = []
  for i in range(len(ix)):
    if (ix[i] < 1):
      ix[i] = 1
    elif (ix[i] > 4095):
      ix[i] = 4095

  for i in range(len(iy)):
    if (iy[i] < 1):
      iy[i] = 1
    elif (iy[i] > 4095):
      iy[i] = 4095

###Create the differential amounts for the correction
  fx = xr - ix
  fy = yr - iy

###Create the arrays to store the corrected values
  xc = np.zeros((len(xr)))
  yc = np.zeros((len(yr)))

###Correct the positions
  for i in range(len(xr)):
    xc[i] = (1-fx[i])*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcxim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcxim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcxim[int(iy[i])][int(ix[i])]

    yc[i] = (1-fx[i])*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])-1] + \
           (1-fx[i])*( fy[i] )*gcyim[int(iy[i])][int(ix[i])-1] + \
           ( fx[i] )*(1-fy[i])*gcyim[int(iy[i])-1][int(ix[i])] + \
           ( fx[i] )*( fy[i] )*gcyim[int(iy[i])][int(ix[i])]

###Return the new arrays
  return xc, yc

################################################################

################################################################
#### Takes in a list of instrumental magnitudes and converts
#### them to normal magnitudes using the exposure time of the
#### observation and a given zeropoint (set in the main script).

def magconvert(magin, imagefile, zp):
  image = fits.open(imagefile)
  exptime = float(image[0].header["EXPTIME"])

  magout = np.zeros((len(magin)))

  for i in range(len(magin)):
    magout[i] = magin[i] + zp + 2.5*np.log10(exptime)

  return magout
  

################################################################

################################################################
#### Transforms input arrays into a format that is used by the 
#### current version of the IDL linear transformation code
#### The "m" indicates it's the master positions and the "p"
#### indicates the positions in the dither to be transformed.

def py2idl(xm,ym,xp,yp,errxm,errxp,errym,erryp):
  lt = len(xm)
  temp = np.zeros(((lt*2),3))
  for i in range(lt):
    temp[i][0] = xp[i]
    temp[i][1] = xm[i]
    temp[i][2] = np.sqrt(errxm[i]**2 + errxp[i]**2)
    temp[i+lt][0] = yp[i]
    temp[i+lt][1] = ym[i]
    temp[i+lt][2] = np.sqrt(errym[i]**2 + erryp[i]**2)

  return temp

################################################################

################################################################
####

def idl2py(file):
  allpos = []
  with open(file,'r') as f:
    for line in f:
      currentline = line.split(' ')
      for i in range(len(currentline)):
        if (currentline[i] == ''):
          continue
        else:
          allpos.append(float(currentline[i]))
  f.close()

  lt = len(allpos)/2
  pos = np.zeros((lt,2))
  for i in range(lt):
    pos[i][0] = allpos[i]
    pos[i][1] = allpos[i+lt]

  return pos
  
################################################################

################################################################
#### Calculates the median deviation for a given input

def mad(array):
  mstd = np.median(abs(array - np.median(array)))
  return mstd

################################################################

################################################################
#### Finds the matches between two given lists of sources within
#### some tolerance and returns a list of the matching indices

def matchlistid(a1, a2, tol, x1, y1, x2, y2):
  id1 = []
  id2 = []
  for k in range(len(a1)):
    matchid = []
    match = 0
    for l in range(len(a2)):
      if (abs(a1[k][x1] - a2[l][x2]) < tol):
        if (abs(a1[k][y1] - a2[l][y2]) < tol):
          match = match+1
          matchid.append(l)

    if (match == 1):
      id1.append(k)
      id2.append(matchid[0])

  return id1, id2

################################################################

################################################################
#### Finds the matches between two given lists of sources within
#### some tolerance and returns a list of the matching indices

def matchlistpos(a1, a2, tol, x1, y1, x2, y2):

  mat1 = []
  mat2 = []

  for k in range(len(a1)):

    diffx = abs(a2[:,x2] - a1[k][x1])
    diffy = abs(a2[:,y2] - a1[k][y1])
    temp2 = a2[(diffx < tol) & (diffy < tol)]

    if (temp2.shape[0] == 1):
      mat1.append(addnumrow(a1[k]))
      mat2.append(addnumrow(temp2[0]))

  return np.asarray(mat1), np.asarray(mat2)

################################################################

################################################################
#### Converts RA in sexagesimal format to degrees

def sex2deg_ra(hr,min,sec):
  deg = 15.0 * (float(hr) + (float(min)/60.0) + (float(sec)/(60.0*60.0)))
  return deg

################################################################

################################################################
#### Converts Dec in sexagesimal format to degrees

def sex2deg_dec(deg,min,sec):
  dec = (abs(float(deg)) + (float(min)/60.0) + (float(sec)/(60.0*60.0)))
  if (deg<0.0):
    dec = dec * (-1.0)

  return dec

################################################################

################################################################
#### Adds a plot to a set of subplots

def addscatter(axarr,x,y,xlab,ylab,mark="+",col="k",uselab=False,lab=""):
  if (uselab):
    axarr.scatter(x,y,marker=mark,color=col,label=lab)
  else:
    axarr.scatter(x,y,marker=mark,color=col,label=lab)
  axarr.set_xlabel(xlab)
  axarr.set_ylabel(ylab)
  axarr.set_aspect('equal')

################################################################

################################################################
#### Converts a given date (year/month/day) into only years

def date2num(date):
  parts = date.split('-')
  num = float(parts[0]) + float(parts[1])/12.0 + float(parts[2])/(12.0*30.0)
  return num

################################################################

################################################################
####

def calcds9vec(array,resscale,imscale,pra,pdec):
  reslen = np.zeros((len(array)))
  resang = np.zeros((len(array)))

  for i in range(len(array)):
    reslen[i] = np.sqrt(array[i][pra]**2 + array[i][pdec]**2) / resscale * imscale
    tempang = np.arctan(abs(array[i][pdec])/abs(array[i][pra]))

    if (array[i][pdec] > 0) & (array[i][pra] > 0):
      vecang = tempang
    elif (array[i][pdec] > 0) & (array[i][pra] < 0):
      vecang = np.pi - tempang
    elif (array[i][pdec] < 0) & (array[i][pra] < 0):
      vecang = np.pi + tempang
    elif (array[i][pdec] < 0) & (array[i][pra] > 0):
     vecang = 2.0*np.pi - tempang

    resang[i] = np.rad2deg(vecang)

  return reslen, resang

################################################################

################################################################
####

def writeds9vec(pos,magn,ang,file,color="red"):
  f = open(file,'w')
  f.write("# Region file format: DS9 version 4.1\n")
  f.write('global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
  f.write("fk5\n")

  for i in range(len(pos)):
    f.write('# vector('+str(pos[i][0])+','+str(pos[i][1])+','+str(magn[i])+'",'+str(ang[i])+') vector=1 width=3\n')

  f.close()

################################################################

################################################################
####

def writeds9errbox(pos,magn,ang,error,file,color="red"):
  f = open(file,'w')
  f.write("# Region file format: DS9 version 4.1\n")
  f.write('global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
  f.write("fk5\n")

  for i in range(len(pos)):

    rapos = pos[i][0] - 3.3*magn[i]*np.cos(np.deg2rad(ang[i]))/(60.0*60.0)
    decpos = pos[i][1] + 1.00*magn[i]*np.sin(np.deg2rad(ang[i]))/(60.0*60.0)

    f.write('box('+str(rapos)+','+str(decpos)+','+str(error[i][0])+'",'+str(error[i][1])+'",'+str(360.0)+') # width=1\n')

  f.close()

################################################################

################################################################
####

def writeds9circ(pos,file,diam=100.0,color="red"):
  f = open(file,'w')
  f.write("# Region file format: DS9 version 4.1\n")
  f.write('global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
  f.write("fk5\n")

  for i in range(len(pos)):
    f.write('circle('+str(pos[i][0])+','+str(pos[i][1])+','+str(diam)+'") # width=3\n')

  f.close()

################################################################

################################################################
####

def pix2wcs_vel(pmx,pmy,angle,scale,baseline,pmra,pmdec):
  pmdec = (-1*np.sin(angle)*pmx - np.cos(angle)*pmy)*scale/baseline
  pmra= -1.0*(np.cos(angle)*pmx - np.sin(angle)*pmy)*scale/baseline


################################################################

################################################################
####

def find_ind(lst, a):
    return [i for i, x in enumerate(lst) if x==a]

################################################################
####

def writeds9point(pos,file,diam=100.0,color="red",type="circle"):
  f = open(file,'w')
  f.write("# Region file format: DS9 version 4.1\n")
  f.write('global color='+color+' dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
  f.write("fk5\n")

  for i in range(len(pos)):
    f.write('point('+str(pos[i][0])+','+str(pos[i][1])+') # point='+type+'\n')

  f.close()

################################################################
####

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

################################################################
####

