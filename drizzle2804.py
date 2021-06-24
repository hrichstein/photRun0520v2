# import numpy as np
# import os
# import shutil
# from drizzlepac import tweakreg
from drizzlepac import astrodrizzle

# from drizzlepac import staticMask, sky, createMedian, ablot, drizCR, adrizzle

# Peg3 F606W
dir = 'PEGASUS-III_f606w/'
f1 = dir + 'jdan46nsq_flc.fits'
f2 = dir + 'jdan46nuq_flc.fits'
f3 = dir + 'jdan46nxq_flc.fits'
f4 = dir + 'jdan46o1q_flc.fits'

# Peg3 F814W
# dir = 'PEGASUS-III_f814w/'
# f1 = dir + 'jdan47m9q_flc.fits'
# f2 = dir + 'jdan47myq_flc.fits'
# f3 = dir + 'jdan47n2q_flc.fits'
# f4 = dir + 'jdan47n6q_flc.fits'

# Psc2 F606W
# dir = 'PISCES-II_f606w/'
# f1 = dir + 'jdan16fkq_flc.fits'
# f2 = dir + 'jdan16fmq_flc.fits'
# f3 = dir + 'jdan16fpq_flc.fits'
# f4 = dir + 'jdan16ftq_flc.fits'

# Psc2 F814W
# dir = 'PISCES-II_f814w/'
# f1 = dir + 'jdan87mxq_flc.fits'
# f2 = dir + 'jdan87mzq_flc.fits'
# f3 = dir + 'jdan87n2q_flc.fits'
# f4 = dir + 'jdan87n6q_flc.fits'

# [f1,f2,f3]
# [f1,f2,f4]
# [f1,f3,f4]
# [f2,f3,f4]

# tweakreg.TweakReg([f1,f2,f3,f4],updatehdr=True)
astrodrizzle.AstroDrizzle([f1,f2,f3,f4],
                          clean=True,
                          configobj='../Hannah_Data/astrodrizzle_new.cfg')

# astrodrizzle.AstroDrizzle([f1,f2,f3],
#                           clean=True,
#                           configobj='../Hannah_Data/astrodrizzle_new.cfg')
#
# astrodrizzle.AstroDrizzle([f1,f2,f4],
#                           clean=True,
#                           configobj='../Hannah_Data/astrodrizzle_new.cfg')
#
# astrodrizzle.AstroDrizzle([f1,f3,f4],
#                           clean=True,
#                           configobj='../Hannah_Data/astrodrizzle_new.cfg')
#
# astrodrizzle.AstroDrizzle([f2,f3,f4],
#                           clean=True,
#                           configobj='../Hannah_Data/astrodrizzle_new.cfg')






# End
