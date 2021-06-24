import numpy as np
import os

upperDir = '/Volumes/Spare Data/Hannah_Data/'

def sex1(image, outfile, expTime,band) :
  #creates a sextractor line e.g sex img.fits -c default.sex -catalog_name something.txt
  dir = '~/seFilesFLCs/'

  if band=='F606W':
      com = [ "sex ", image, " -c defaultVBand_ACS.sex -GAIN ", str(1./expTime),
      " -CATALOG_NAME ",dir+outfile]

  elif band=='F814W':
      com = [ "sex ", image, " -c defaultIBand_ACS.sex -GAIN ", str(1./expTime),
       " -CATALOG_NAME ",dir+outfile]

  s0=''
  com = s0.join(com)

  return com


def makeSex(nameslist,band='F606W'):

    img_arr = np.genfromtxt(nameslist,usecols=(0,1,4),dtype=str,comments='#')

    if band=='F606W':
        bandir = 'f606w'
    elif band=='F814W':
        bandir = 'f814w'

    for ii in range(len(img_arr[:,0])):

        img_name = img_arr[:,0][ii]
        obj_name = img_arr[:,1][ii]

        # filt1 = img_arr[:,2][ii]
        # filt2 = img_arr[:,3][ii]

        filename = band + 'comsList_' + obj_name + '.txt'
        dir = obj_name + '_' + bandir + '/'

        if ii==0:
            com_file = open(dir+filename,'w')
        elif obj_name==img_arr[:,1][ii-1]:
            com_file = open(dir+filename,'a')
        else:
            com_file = open(dir+filename,'w')

        img_tag = img_name.replace('_WJ2.fits','')
        outfile = 'se_' + img_tag + '_' + str(obj_name) + "_" + band + '.dat'

        expTime_str = img_arr[:,2][ii]
        expTime_fl = np.float(expTime_str)

        tmpCom = sex1(img_name,outfile,expTime_fl,band)

        com_file.write('%s\n' % tmpCom)

    com_file.close()

    return None
