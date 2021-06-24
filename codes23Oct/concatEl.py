import numpy as np

catDir = 'drcPhot10Nov/'


tri_arr = ['TRIANGULUM-II-EAST','TRIANGULUM-II-WEST']
tuc_arr = ['TUCANA-II-NE','TUCANA-II-SE']
urs_arr = ['URSA-MAJOR-II-WEST','URSA-MAJOR-II-EAST']
tu4_arr = ['TUCANA-IV-SOUTH','TUCANA-IV-NORTH']
tu3_arr = ['TUCANA-III-WEST','TUCANA-III-EAST']
seg_arr = ['SEGUE-1-WEST','SEGUE-1-EAST']
hor_arr = ['HOROLOGIUM-II-WEST','HOROLOGIUM-II-EAST']
boo_arr = ['BOOTES-II-SOUTH','BOOTES-II-NORTH']


def concat(name_arr,source,catDir='../forElena16Nov/forElena31May/'):

    for ff in range(len(name_arr)):

        if ff==0:
            f_out = np.genfromtxt(catDir + name_arr[ff] + '_cat31May.dat')
            fhead = np.genfromtxt(catDir + name_arr[ff] + '_cat31May.dat',
                                  names=True)
        else:
            f_temp = np.genfromtxt(catDir+ name_arr[ff] + '_cat31May.dat')
            f_out = np.vstack((f_out,f_temp))

    colNs = np.array(fhead.dtype.names)
    header = ' '.join(colNs)

    np.savetxt(catDir + source + '_cat31May.dat',f_out,header=header,
               fmt='%1.5f')

    return None

#
# concat(tri_arr,'TRIANGULUM-II')
# concat(tuc_arr,'TUCANA-II')
concat(urs_arr,'URSA-MAJOR-II')
concat(tu4_arr,'TUCANA-IV')
concat(tu3_arr,'TUCANA-III')
concat(seg_arr,'SEGUE-1')
concat(hor_arr,'HOROLOGIUM-II')
concat(boo_arr,'BOOTES-II')
