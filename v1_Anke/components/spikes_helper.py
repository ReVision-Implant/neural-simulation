import pandas as pd
from hdf5 import HDF5
import matplotlib.pyplot as plt
import numpy as np
import os

def get_spikes(nodes_dirs, spikes_dirs, spikes_bg_dirs, radius=None, v1=True, **kwargs):
    
    nodes_dir = nodes_dirs[0] if isinstance(nodes_dirs, list) else nodes_dirs
    spikes_dir = spikes_dirs[0] if isinstance(spikes_dirs, list) else spikes_dirs 

    if v1 == True:
        node_pos = HDF5(nodes_dir).get_positions_v1()
    elif v1 == False:
        node_pos = HDF5(nodes_dir).get_positions()

    n_spikes = np.zeros(np.shape(node_pos)[0])

    spikes = pd.read_csv(spikes_dir, sep='\s+')
    for ind in spikes.index:
        if spikes['timestamps'][ind] < 100:
            n_spikes[spikes['node_ids'][ind]] += 1

    if spikes_bg_dirs is not None:
        spikes_bg_dir = spikes_bg_dirs[0] if isinstance(spikes_bg_dirs, list) else spikes_bg_dirs 
        spikes_bg = pd.read_csv(spikes_bg_dir, sep='\s+')
        for ind in spikes_bg.index:
            if spikes_bg['timestamps'][ind] < 100:
                n_spikes[spikes_bg['node_ids'][ind]] = max(0, n_spikes[spikes_bg['node_ids'][ind]] - 1)
    
    if radius is not None:
        circle = node_pos[:,0]**2 + node_pos[:,2]**2 < radius**2

    return node_pos[circle,:], n_spikes[circle]

def get_path_kwargs(exp, patterns, amplitudes, networks):

    # Initialisation
    kwargs = dict(nodes_dirs=[], spikes_dirs=[], spikes_bg_dirs=[], electrodes_dirs=[])
    path = os.path.dirname(os.path.abspath(__file__))
    path_to_v1 = os.path.dirname(path)

    # Convert strings to lists if necessary
    patterns = [patterns] if not isinstance(patterns, list) else patterns
    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    networks = [networks] if not isinstance(networks, list) else networks

    # Convert integers to strings if necessary
    exp = 'exp' + str(exp) if str(exp)[:3] != 'exp' else str(exp)
    patterns = ['pattern'+str(i) for i in patterns]
    amplitudes = [str(i)+'uA' for i in amplitudes]
    networks = ['network'+str(i) for i in networks]

    # Iterate over all parameters and get paths for each combination
    for amplitude in amplitudes:
        for pattern in patterns:
            for network in networks:

                nodes_dir = os.path.join(path_to_v1,'networks', network, 'v1_nodes.h5')
                spikes_dir = os.path.join(path_to_v1, exp, 'output', pattern, amplitude, network, 'spikes.csv')
                spikes_bg_dir = os.path.join(path_to_v1, exp, 'output', 'bkg', network, 'spikes.csv')
                electrodes_dir = os.path.join(path, 'electrodes', exp, pattern, 'electrodes.csv')

                kwargs['nodes_dirs'].append(nodes_dir)
                kwargs['spikes_dirs'].append(spikes_dir)
                kwargs['spikes_bg_dirs'].append(spikes_bg_dir)
                kwargs['electrodes_dirs'].append(electrodes_dir)

    return kwargs

kwargs = get_path_kwargs(1,3,10,1)
print(kwargs)
# spikes = get_spikes(**kwargs)
# print(spikes)

def print_args(names=False, **kwargs):
    print(kwargs)
    length = 1
    for key, value in kwargs.items():
        length *= len(value)
    res = ['']*length
    period = length
    for key, value in kwargs.items():
        period /= len(value)
        iters = length/period
        for i in range(length):
            value_i = (i // period) % len(value)
            value_i = int(value_i)
            if names:
                res[i] += key + str(value[value_i]) + '_'
            else:
                res[i] += str(value[value_i]) + '_'
    for i in range(length):
        res[i] = res[i][:-1]

    print(res)

kwargs = {'electrodes':[1], 'stim_type':'-', 'amplitudes':[10,20,30], 'networks':[0,1]}
print_args(**kwargs, names=False)