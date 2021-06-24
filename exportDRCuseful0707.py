import numpy as np

drc_dir  = '/Users/hr8jz/Box Sync/Research/source_lists/june13/'
# workDir='./catRawMags1305/catDir/'

def exp_DRC(targname):

    workDir = './catRawMags1305/catDir_'+targname+'/'

    drc = np.genfromtxt(drc_dir+targname+'_sfErr.dat',names=True)

    outArr = np.zeros((len(drc),13))

    outArr[:,0] = drc['RA_v']
    outArr[:,1] = drc['DEC_v']
    outArr[:,2] = drc['x_v']
    outArr[:,3] = drc['y_v']
    outArr[:,4] = drc['magRaw_v']
    outArr[:,5] = drc['class_Star_v']
    outArr[:,6] = drc['RA_i']
    outArr[:,7] = drc['DEC_i']
    outArr[:,8] = drc['x_i']
    outArr[:,9] = drc['y_i']
    outArr[:,10] = drc['magRaw_i']
    outArr[:,11] = drc['class_Star_i']
    outArr[:,12] = np.arange(0,len(drc),1)

    header = 'RA_f606w DEC_f606w x_f606w y_f606w magr_f606w c_star_f606w RA_f814w DEC_f814w x_f814w y_f814w magr_f814w c_star_f814w id'
    form = "%1.7f %1.7f %1.5f %1.5f %1.5f %1.3f %1.7f %1.7f %1.5f %1.5f %1.5f %1.3f %d"

    np.savetxt(workDir+'drc_useful_'+targname+'.dat',outArr,header=header,fmt=form)

    return None

# targname_arr = ['HYDRA-II','PEGASUS-III','PHOENIX-II','RETICULUM-II','TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE',
# 'TUCANA-II-NW','TUCANA-II-SE','TUCANA-II-SW','SAGITTARIUS-II']
# targname_arr=['HOROLOGIUM-I']
# for c1,tt in enumerate(targname_arr):
#     exp_DRC(tt)


# end
