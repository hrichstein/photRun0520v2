import numpy as np
import matplotlib.pyplot as plt

def makeCMD(targname,dir='./',filt='F606W'):

    match = np.genfromtxt(dir+'newOLDmatch_FULL_'+filt+'.dat',names=True)

    idx = np.logical_and(np.logical_and(match['std_f606wO']<=0.5, match['std_f814wO']<=0.5),np.logical_and(match['std_f606wN']<=0.5, match['std_f814wN']<=0.5))

    match = match[idx]

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,6.5),sharex=True,sharey=True)

    ax1.scatter(match['magr_f606wO']-match['magr_f814wO'],match['magr_f606wO'],s=5,label='OLD')
    ax1.scatter(match['magr_f606wN']-match['magr_f814wN'],match['magr_f606wN'],s=5,label='NEW')

    ax2.scatter(match['magr_f606wO']-match['magr_f814wO'],match['magr_f814wO'],s=5,label='OLD')
    ax2.scatter(match['magr_f606wN']-match['magr_f814wN'],match['magr_f814wN'],s=5,label='NEW')

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

    plt.savefig(dir+targname+'_OldNew.png',dpi=600,bbox_inches='tight')

    plt.close()

    return None
