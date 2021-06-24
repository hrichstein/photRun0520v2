import numpy as np
import matplotlib.pyplot as plt

def makeCMD(targname,dir='./'):

    match = np.genfromtxt(dir+targname+'_APER2flcMatch_1408.dat',names=True)

    idx = np.logical_and(np.logical_and(match['stdev_f606w']<=0.5, match['stdev_f814w']<=0.5),np.logical_and(match['s606']<=0.5, match['s814']<=0.5))

    match = match[idx]

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)

    ax1.scatter(match['m606c']-match['m814c'],match['m606c'],s=5,label='APER')
    ax1.scatter(match['magZPT_f606w']-match['magZPT_f814w'],match['magZPT_f606w'],s=5,label='FLC')

    ax2.scatter(match['m606c']-match['m814c'],match['m814c'],s=5,label='APER')
    ax2.scatter(match['magZPT_f606w']-match['magZPT_f814w'],match['magZPT_f814w'],s=5,label='FLC')

    ax1.set_ylim(28,17.5)
    ax1.set_xlim(-1.5,1.5)
    ax1.set_xlabel('F606W-F814W')
    ax1.set_ylabel('F606W')

    ax2.set_xlabel('F606W-F814W')
    ax2.set_ylabel('F814W')

    ax1.set_title(targname)
    ax2.set_title(targname)

    ax1.legend()
    ax2.legend()

    plt.savefig(dir+targname+'_APERflc.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None
