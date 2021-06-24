import numpy as np
import matplotlib.pyplot as plt

def f2mag_dirs(targname,date='20Aug',workDir='./'):

    return workDir+'catRawMags'+date+'/catDir_'+targname+'/'

def outSGcut(namefile,sg_cut=0.5):

    file = open('catRawMags20Aug/sgCut_everyWerr.dat','w')
    # file = open('catRawMags20Aug/sgCut_MWwErr.dat','w')
    file.write('# f606w f814wn errf606w errf814w \n')
    targname_arr = np.genfromtxt(namefile,dtype='str')
    # targname_arr = ['HOROLOGIUM-I','PHOENIX-II','RETICULUM-II']
    # targname_arr = ['HYDRA-II','PEGASUS-III','SAGITTARIUS-II',\
    # 'TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE',\
    # 'TUCANA-II-SE', 'HOROLOGIUM-I','PHOENIX-II','RETICULUM-II']

    # fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)
    for tt,tname in enumerate(targname_arr):

        dir = f2mag_dirs(tname,date='20Aug',workDir='./')
        # flc606 = np.genfromtxt(dir+'classStarCat_F606Wmatch_'+tname+'.dat',names=True)
        #
        # flc814 = np.genfromtxt(dir+'classStarCat_F814Wmatch_'+tname+'.dat',names=True)
        #
        # len606 = len(flc606)
        # len814 = len(flc814)
        #
        # if len606>len814:
        #     flc = flc606
        #     print(tname,'F606W')
        # else:
        #     flc = flc814
        #     print(tname,'F814W')
        # flc = np.genfromtxt(dir+'sgCut_'+tname+'_tcCut.dat',names=True)
        flc = np.genfromtxt(dir+'sgCut_'+tname+'.dat',names=True)

        colNs = np.array(flc.dtype.names)

        flc_idx = np.logical_and(np.logical_and(flc['magZPT_f606w']>=19.5,flc['c_star_f606w']>=sg_cut),flc['c_star_f814w']>=sg_cut)

        flc_g = flc[flc_idx]

        for ff, line in enumerate(flc_g):
            file.write('{0} {1} {2} {3} \n'.format(flc_g['magAbs_f606w'][ff],flc_g['magAbs_f814w'][ff],
            flc_g["magZPTerr_f606w"][ff],flc_g["magZPTerr_f814w"][ff]))


        # s0 = ' '
        # header = s0.join(colNs)
    # header = 'f606w f814w'
    # flc_out = np.array([f606w,f814w]).flatten()
    # print(type(flc_out))
        # np.savetxt(dir+'sgCut_'+tname+'.dat',flc_g,header=header)
    file.close()
    # np.savetxt(dir+'sgCut_ALL.dat',flc_out.T,header=header)

    return None

def main():

    outSGcut('targnamesDirections2.txt',sg_cut=0.95)
    # outSGcut('targnamesPost.txt',sg_cut=0.95)


if __name__ == '__main__':
    main()
