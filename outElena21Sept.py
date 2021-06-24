import numpy as np
import matplotlib.pyplot as plt

def f2mag_dirs(targname,date='20Aug',workDir='./'):

    return workDir+'catRawMags'+date+'/catDir_'+targname+'/'

def outSGcut(namefile,sg_cut=0.5):

    # file = open('catRawMags20Aug/sgCut_all.dat','w')
    # file = open('catRawMags20Aug/sgCut_LMC.dat','w')
    # file.write('# f606w f814wn\n')
    targname_arr = np.genfromtxt(namefile,dtype='str')
    # targname_arr = ['HOROLOGIUM-I','PHOENIX-II','RETICULUM-II']
    # targname_arr = ['HYDRA-II','PEGASUS-III','SAGITTARIUS-II',\
    # 'TRIANGULUM-II-EAST','TRIANGULUM-II-WEST','TUCANA-II-NE',\
    # 'TUCANA-II-SE']

    # fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)
    for tt,tname in enumerate(targname_arr):

        dir = f2mag_dirs(tname,date='20Aug',workDir='./')
        flc606 = np.genfromtxt(dir+'classStarCat_F606Wmatch_'+tname+'.dat',names=True)

        flc606c = np.genfromtxt(dir+'classStarCat_F606Wmatch_'+tname+'.dat')

        flc814 = np.genfromtxt(dir+'classStarCat_F814Wmatch_'+tname+'.dat',names=True)

        flc814c = np.genfromtxt(dir+'classStarCat_F814Wmatch_'+tname+'.dat')

        len606 = len(flc606)
        len814 = len(flc814)

        if len606>len814:
            flc = flc606c
            flcd = flc606
            print(tname,'F606W')
        else:
            flc = flc814c
            flcd = flc814
            print(tname,'F814W')
        # flc = np.genfromtxt(dir+'sgCut_'+tname+'_tcCut.dat',names=True)

        colNs = np.array(flcd.dtype.names)

        xt_606 = np.int(np.where(colNs=='xt1_f606w')[0])
        yt_606 = np.int(np.where(colNs=='yt1_f606w')[0])

        xt_814 = np.int(np.where(colNs=='xt1_f814w')[0])
        yt_814 = np.int(np.where(colNs=='yt1_f814w')[0])

        mean_f606w = np.int(np.where(colNs=='mean_f606w')[0])
        mean_f814w = np.int(np.where(colNs=='mean_f814w')[0])

        stdev_f606w = np.int(np.where(colNs=='stdev_f606w')[0])
        stdev_f814w = np.int(np.where(colNs=='stdev_f814w')[0])

        magZPT_f606w = np.int(np.where(colNs=='magZPT_f606w')[0])
        magZPT_f814w = np.int(np.where(colNs=='magZPT_f814w')[0])

        magZPTerr_f606w = np.int(np.where(colNs=='magZPTerr_f606w')[0])
        magZPTerr_f814w = np.int(np.where(colNs=='magZPTerr_f814w')[0])

        c_star_f606w = np.int(np.where(colNs=='c_star_f606w')[0])
        c_star_f814w = np.int(np.where(colNs=='c_star_f814w')[0])

        flc_idx = np.logical_and(flcd['c_star_f606w']>=sg_cut,flcd['c_star_f814w']>=sg_cut)

        flc_g = flc[flc_idx]

        # outArr = np.zeros((len(flc_g),14))
        #
        # outArr[:,0] = flc_g[:,xt_606]

        outArr = flc_g[:,[xt_606,yt_606,mean_f606w,stdev_f606w,magZPT_f606w,\
        magZPTerr_f606w,c_star_f606w,xt_814,yt_814,mean_f814w,stdev_f814w,\
        magZPT_f814w,magZPTerr_f814w,c_star_f814w]]

        # outArr = flc_g[:,[xt_606,yt_606]]

        header = 'xt_f606w yt_f606w mean_f606w stdev_f606w magZPT_f606w '
        header += 'magZPTerr_f606w c_star_f606w xt_f814w yt_f814w mean_f814w '
        header += 'stdev_f814w magZPT_f814w magZPTerr_f814w c_star_f814w'

        np.savetxt('catsSept21/catsSept21_'+tname+'.dat',outArr,header=header)
        # for ff, line in enumerate(flc_g):
        #     file.write('{0} {1} \n'.format(flc_g['magAbs_f606w'][ff],flc_g['magAbs_f814w'][ff]))


        # s0 = ' '
        # header = s0.join(colNs)
    # header = 'f606w f814w'
    # flc_out = np.array([f606w,f814w]).flatten()
    # print(type(flc_out))
        # np.savetxt(dir+'sgCut_'+tname+'.dat',flc_g,header=header)
    # file.close()
    # np.savetxt(dir+'sgCut_ALL.dat',flc_out.T,header=header)

    return None

def main():

    # makeNewCMD('targnamesDirections2.txt',sg_cut=0.95)
    outSGcut('targnamesPost.txt',sg_cut=0.5)


if __name__ == '__main__':
    main()
