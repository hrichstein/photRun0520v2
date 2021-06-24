import numpy as np

def runSEcoms(targnames):

    filters = ['F606W','F814W']
    lil_f = ['f606w','f814w']

    names = np.genfromtxt(targnames,dtype=str)

    com_file = open('masterList.sh','w')

    for nn in range(len(names)):
        for ff in range(len(filters)):
            com1 = 'cp mountFolder/' + names[nn] + '_' + lil_f[ff] + '/' + filters[ff] + 'comsList_' + names[nn] + '.txt . '

            com_file.write('%s\n' % com1)

            com2 = 'cp mountFolder/' + names[nn] + '_' + lil_f[ff] + '/crClean/*WJ2.fits . '

            com_file.write('%s\n' % com2)

            com3 = 'python run_sex.py ' + filters[ff] + 'comsList_' + names[nn] + '.txt '

            com_file.write('%s\n' % com3)

            com4 = 'rm *WJ2.fits '

            com_file.write('%s\n' % com4)

            com5 = 'rm ' + filters[ff] + 'comsList_' + names[nn] + '.txt '
            com_file.write('%s\n' % com5)

    com_file.close()

    return None
