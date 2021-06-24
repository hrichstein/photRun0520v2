dir = './'

def matchWJCs(targname,filt,dir='./',matchtol=5):

    # xt, yt = 11,12
    # magr,id = 7,8

    jdanUse = getJdan(targname,filt)
    outName = "master_ids_"+targname+"_"+filt+".dat"

    suffix = '_oc.dat'

    master = np.genfromtxt(dir+jdanUse[0]+"_"+targname+'_'+filt +suffix,names=True)
    masterCat = np.loadtxt(dir+jdanUse[0]+"_"+targname+'_'+filt +suffix)

    colNs = np.array(master.dtype.names)


    xt = np.int(np.where(colNs=='xo')[0])
    yt = np.int(np.where(colNs=='yo')[0])
    xtstr = 'xo'
    ytstr = 'yo'

    magr = np.int(np.where(colNs=='magr')[0])
    id = np.int(np.where(colNs=='id')[0])
    # Create an array of zeros with columns equal to the number of non-master dithers to store the matching id for each source
    matchids = np.zeros((len(master), (len(jdanUse)-1)))
    master = np.hstack((masterCat, matchids))

    # Loop through other images
    for dd in range(len(jdanUse)-1):
        # Load catalogs
        cat = np.genfromtxt(dir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix,names=True)
        catCat = np.loadtxt(dir+jdanUse[dd+1]+"_"+targname+'_'+filt+suffix)

        colCNs = np.array(cat.dtype.names)

        nF = True
        row = 0

        while (nF): # not finished
            matchrows = cat[(abs(master[row][xt] - cat[xtstr]) <= matchtol) & (abs(master[row][yt] - cat[ytstr]) <= matchtol)]

    #         # Setting the proper column number to the matching index.
            if (len(matchrows) == 1):
              master[row][xt+dd+2] = matchrows[0][id]
              row = row + 1

            elif (len(matchrows) > 1):
                magDif = np.zeros((len(matchrows),1))
                for mm in range(len(matchrows)):
                    magDif[mm] = master[row][magr] - matchrows[mm][magr]
                small = np.argmin(magDif)
                master[row][xt+dd+2] = matchrows[small][id]
                row += 1

            else:
              master = np.delete(master, row, 0)

            if (row >= len(master)):
                nF = False

    header =  "flags RA DEC xr yr flux c_star magr id xc yc xt yt id2 id3 id4"
    form = "%d %1.7f %1.7f %1.4f %1.4f %1.4f %1.3f %1.4f %d %1.4f %1.4f %1.4f %1.4f %d %d %d"

    np.savetxt(dir+outName,master, header=header, fmt=form)


    return None

########################
# Not a function right here - quick and dirty code for HOR-I F814W

d1_cat = np.genfromtxt(dir+'jdan21l8q_HOROLOGIUM-I_F814W_oc.dat')
d2_cat = np.genfromtxt(dir+'jdan21laq_HOROLOGIUM-I_F814W_oc.dat')
d3_cat = np.genfromtxt(dir+'jdan21lhq_HOROLOGIUM-I_F814W_oc.dat')
d4_cat = np.genfromtxt(dir+'jdan21llq_HOROLOGIUM-I_F814W_oc.dat')
cat_arrA = np.array([d1_cat,d2_cat,d3_cat,d4_cat])

# Replace master_ids_HOROLOGIUM-I_F814W.dat with the output file from matchWJCs
# if the name differs

master = np.genfromtxt(dir+"master_ids_HOROLOGIUM-I_F814W.dat",names=True)
masterCat = np.loadtxt(dir+"master_ids_HOROLOGIUM-I_F814W.dat")

colNs = np.array(master.dtype.names)

ra_id = np.int(np.where(colNs=='RA')[0])
dec_id = np.int(np.where(colNs=='DEC')[0])
flu_id = np.int(np.where(colNs=='flux')[0])
fla_id = np.int(np.where(colNs=='flags')[0])
cs_id = np.int(np.where(colNs=='c_star')[0])
magr = np.int(np.where(colNs=='magr')[0])

xr = np.int(np.where(colNs=='xr')[0])
yr = np.int(np.where(colNs=='yr')[0])
xc = np.int(np.where(colNs=='xc')[0])
yc = np.int(np.where(colNs=='yc')[0])

id2 = np.int(np.where(colNs=='id2')[0])
id3 = np.int(np.where(colNs=='id3')[0])
id4 = np.int(np.where(colNs=='id4')[0])
# xt, yt are xo,yo in dithers 2-4 as the transformed positions
# don't change for the first dither
xt = np.int(np.where(colNs=='xt')[0])
yt = np.int(np.where(colNs=='yt')[0])

id = np.int(np.where(colNs=='id')[0])
coordRows = masterCat[:,[ra_id,dec_id,flu_id,fla_id,cs_id]]

nCo = 4*int(9) # 4 is number of dithers (normally have a different variable, not just the number)
newCols = np.zeros((len(coordRows), nCo))

# rowsMast = np.transpose(masterCat)

jj = 0
cc = 0
while jj < 4:
    cat = cat_arrA[jj]

    if jj==0:
        idcol = id
    elif jj==1:
        idcol = id2
    elif jj==2:
        idcol = id3
    elif jj==3:
        idcol = id4

    newIDcol = masterCat[:,idcol]
    idx = np.asarray(newIDcol,int)

    reg = cat[idx]

    newCols[:,cc] = reg[:,magr]
    newCols[:,cc+jj+4] = reg[:,ra_id]
    newCols[:,cc+jj+5] = reg[:,dec_id]

    newCols[:,cc+jj+12] = reg[:,xr]
    newCols[:,cc+jj+13] = reg[:,yr]

    newCols[:,cc+jj+20] = reg[:,xc]
    newCols[:,cc+jj+21] = reg[:,yc]

    newCols[:,cc+jj+28] = reg[:,xt]
    newCols[:,cc+jj+29] = reg[:,yt]

    cc += 1
    jj += 1


magList = np.hstack((coordRows, newCols))

header = 'RA DEC flux flags c_star mag1 mag2 mag3 mag4 ra1 dec1 ra2 dec2 ra3 dec3 ra4 dec4 xr1 yr1 xc1 yc1 xt1 yt1 xr2 yr2 xc2 yc2 xt2 yt2 xr3 yr3 xc3 yc3 xt3 yt3 xr4 yr4 xc4 xc4 xt4 yt4'
form = '%1.7f %1.7f %1.4f %d %1.3f %1.4f %1.4f %1.4f %1.4f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.7f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f'

np.savetxt(dir+'matched_w_MagsPos.dat',magList,header=header,fmt=form)
