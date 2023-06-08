import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
import pandas as pd
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../../bio_components')
sys.path.append('../bio_components')
from bio_components.helper import get_params, get_spikes

# Initialisation
# assert isinstance(electrodes, list) + isinstance(amplitudes, list) + isinstance(exp, list) == 1 # Correlation between either different electrodes or different amplitudes
names = []
n_spikes = []

exp = [0,3]
electrodes = [1]
amplitudes = [10]
networks = [0,1,2]
stim_type = 'g'

corrs = np.zeros((len(exp)*len(electrodes)*len(amplitudes), len(exp)*len(electrodes)*len(amplitudes)))
N =  np.zeros((len(exp)*len(electrodes)*len(amplitudes), len(exp)*len(electrodes)*len(amplitudes)))

# Iterate over electrodes/amplitudes: get spikes for all networks and get corrcoef with spike lists of previous iterations
for i, exp in enumerate(exp):
    for j,electrode in enumerate(electrodes):
        for k,amplitude in enumerate(amplitudes):
            if stim_type in [0,'0']:
                exp = 0
                stim_type = 'g'
                electrode_temp = 0
            else:
                electrode_temp = electrode
                exp = 3
            n_spikes_temp = np.array([])
            for network in networks:
                n_spikes_temp = np.append(n_spikes_temp, get_spikes(**get_params(exp, electrode, stim_type, amplitude, network))[1])
            n_spikes.append(n_spikes_temp)
            n = i*len(electrodes)*len(amplitudes)+j*len(amplitudes)+k
            for u in range(n):
                corrs[u,n] = np.corrcoef(n_spikes[n],n_spikes[u])[0,1]
                N[u,n] = np.sum((n_spikes[u]+n_spikes[n])>0)
                names.append(str(n)+'_'+str(u))

print(corrs)
