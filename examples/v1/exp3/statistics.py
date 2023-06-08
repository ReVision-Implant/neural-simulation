import pandas as pd 
import sys
import numpy as np
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import get_corr_inter, get_corr_intra, get_p_values
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def corr_anal(sigma, amplitude):
    exp = '3'
    stim_type = '-'
    electrodes = ['1','2','3']
    networks = ['0','1','2']
    intra_corrs = []
    inter_corrs = []


    for i, electrode in enumerate(electrodes):
        corrs, _ = get_corr_intra(sigma, exp, electrode, stim_type, amplitude, networks)
        intra_corrs.append(np.mean(corrs)) # average over pairwise correlations
        for electrode2 in electrodes[:i]:
            corrs, _ = get_corr_inter(exp, [electrode, electrode2], stim_type, amplitude, networks)
            inter_corrs.append(corrs)
    intra_corrs = np.mean(intra_corrs) # average over electrodes

    return intra_corrs, np.array(inter_corrs).flatten()

def sigma(): 
    intra_corrs = np.zeros(20)
    inter_corrs = np.zeros((20,3))
    sigmas = range(0,20)
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True)
    for i, amplitude in enumerate(range(10,40,10)):
        for sigma in sigmas:
            intra_corr, inter_corr = corr_anal(sigma, amplitude)
            print('sigma', sigma, ':', intra_corr, inter_corr)
            intra_corrs[sigma] = intra_corr
            inter_corrs[sigma,:] = inter_corr

        axs[i].plot(sigmas, intra_corrs)
        axs[i].plot(sigmas, inter_corrs[:,0])
        axs[i].plot(sigmas, inter_corrs[:,1])
        axs[i].plot(sigmas, inter_corrs[:,2])
        axs[i].set_title('I=' +str(amplitude) + 'uA')
        axs[i].set_ylabel('correlation')
        axs[i].grid(axis='y', color='lightgrey')
    axs[2].set_xlabel('$\sigma$ of gaussian filter')
    axs[0].legend(['intra', 'inter12', 'inter13', 'inter23'])
    plt.tight_layout()
    plt.savefig('exp3/sigma.png')
    plt.show()

if __name__ == "__main__":

    inter_corr, _, N = get_corr_inter(['3','3-'],1,'-',30,[0,1,2])
    # intra_corr, _ = get_corr_intra(10,['3','3-'],1,'-',30,[0,1,2])
    # intra_corr = np.mean(intra_corr)
    print(N)
    p_values = get_p_values(inter_corr, 0.90, N)
    print(p_values)
