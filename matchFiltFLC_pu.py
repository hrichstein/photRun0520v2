import numpy as np
# import matplotlib.pyplot as plt

def matchFilt(targname,dir='./',matchtol=3,newdir='./'):

    f606wN = np.genfromtxt(dir+'magZPTedAll_F606W.dat',names=True)
    f606w = np.genfromtxt(dir+'magZPTedAll_F606W.dat')

    # f814wN = np.genfromtxt(dir+'magZPTedAll_F814W.dat',names=True)
    # f814w = np.genfromtxt(dir+'magZPTedAll_F814W.dat')
    f814wN = np.genfromtxt(newdir+targname + "_filtTrans_pu.dat",names=True)
    f814w = np.genfromtxt(newdir+targname + "_filtTrans_pu.dat")

    id606 = np.zeros((len(f606w),1))
    id606[:,0] = np.arange(0,len(f606w),1)

    f606w = np.hstack((f606w,id606))
    col606 = np.array(f606wN.dtype.names)

    id814 = np.zeros((len(f814w),1))
    id814[:,0] = np.arange(0,len(f814w),1)

    f814w = np.hstack((f814w,id814))
    col814 = np.array(f814wN.dtype.names)

    # Getting columns of x,y of transformed F814W to F606W
    xI = np.int(np.where(col814=='x_f606wTrans')[0])
    yI = np.int(np.where(col814=='y_f606wTrans')[0])
    magI = np.int(np.where(col814=='magZPT_wa')[0])
    idI = len(col814)

    xV = np.int(np.where(col606=='xt1')[0])
    yV = np.int(np.where(col606=='yt1')[0])
    magV = np.int(np.where(col606=='magZPT_wa')[0])
    idV = len(col606)

    master_in = f606w[:,[xV,yV,magV,idV]]
    x,y,magrV,idV_mas = 0,1,2,3

    cat = f814w
    matchids = np.zeros((len(master_in),1))

    nF_out = True
    matchtol=matchtol

    len606 = len(f606w)
    len814 = len(f814w)

    minLen = np.min([len606,len814])

    while nF_out:

        master, matchids = matchlistID(master_in,cat,matchtol,x,y,xI,yI,idI)

        if len(master)>=int(0.75*minLen):
            nF_out = False
            print('Minimum Number Reached: %d' % len(master),targname)

        else:
            print('Need More Stars')
            print("Pixel Tolerance: %d, Number Stars: %d" % (matchtol,len(master)))
            matchtol += 1
            if matchtol <= 10:
                master_in = f606w[:,[xV,yV,magV,idV]]
                matchids = np.zeros((len(master_in),1))
            else:
                print("Sacrificing number of stars for quality of matches.")
                nF_out = False

    master = np.hstack((master,matchids))
    print(targname, len(master)/minLen)

    xV_mas, yV_mas, magV_mas, idV_mas, idI_mas = 0, 1, 2, 3, 4

    # newCols = np.zeros((len(master),3))

    idCol606 = master[:,idV_mas]
    idx606 = np.asarray(idCol606,int)
    reg606 = f606w[idx606]

    idCol814 = master[:,idI_mas]
    idx814 = np.asarray(idCol814,int)
    reg814 = f814w[idx814]

    outArr = np.hstack((reg606,reg814))

    col606_temp = np.array(["{}{}".format(n,'_f606w') for n in col606])
    s0 = ' '
    header606 = s0.join(col606_temp)
    header606 += ' id_f606w '

    col814_temp = np.array(["{}{}".format(n,'_f814w') for n in col814])
    header814 = s0.join(col814_temp)
    header814 += ' id_f814w'

    header = header606 + header814

    np.savetxt(newdir+targname+'_allMatchedZPTed_pu.dat',outArr,header=header)
    print(newdir+targname+'_allMatchedZPTed_pu.dat')

    return None

def matchlistID(master,cat,matchtol,x1,y1,x2,y2,id_mat):

    matchids_in = np.zeros((len(master),1))

    nF = True
    row = 0

    while nF:

        matchrows = cat[(abs(master[row][x1] - cat[:,x2]) \
            <= matchtol) & (abs(master[row][y1] - cat[:,y2])<= matchtol)]

        if (len(matchrows) == 1):
            matchids_in[row][0] = matchrows[0][id_mat]
            row += 1

        elif (len(matchrows) > 1):
            distDiff = np.zeros((len(matchrows),1))
            for dd in range(len(matchrows)):
                distDiff[dd] = np.sqrt( (master[row][x1] - \
                matchrows[dd][x2])**2 +  (master[row][y1] \
                                    - matchrows[dd][y2])**2)
            small = np.argmin(distDiff)
            matchids_in[row][0] = matchrows[small][id_mat]
            row += 1

        else:
            master = np.delete(master,row,0)
            matchids_in = np.delete(matchids_in,row,0)

        if (row >= len(master)):
            u, udx = np.unique(matchids_in,return_index=True)

            if len(udx)<len(master):
                master = master[udx]
                matchids_in = matchids_in[udx]
                nF = False

            elif len(udx)==len(master):
                nF = False


    return master,matchids_in

#
