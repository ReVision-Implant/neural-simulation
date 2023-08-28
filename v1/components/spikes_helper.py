import pandas as pd
from hdf5 import HDF5
import matplotlib.pyplot as plt
import numpy as np

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


def get_params(exp, electrodes, stim_type, amplitudes, networks):

    # Initialisation
    nodes = []
    spikes = []
    spikes_bg = []
    elec_loc = []

    # Convert strings to lists if necessary
    electrodes = [electrodes] if not isinstance(electrodes, list) else electrodes
    amplitudes = [amplitudes] if not isinstance(amplitudes, list) else amplitudes
    networks = [networks] if not isinstance(networks, list) else networks
       

    # Convert integers to strings if necessary
    exp = str(exp)
    electrodes = [str(i) for i in electrodes]
    amplitudes = [str(i) for i in amplitudes]
    networks = [str(i) for i in networks]

    # Networks to use
    density = '25' if exp=='1' else '100'
    radius = 400 if exp=='1' else 200

    # Iterate over all parameters and get paths for each combination
    for amplitude in amplitudes:
        for electrode in electrodes:
            for network in networks:
                nodes.append('networks_' + density + '/network' + network + '/v1_nodes.h5')
                spikes.append('exp' + exp + '/output/'+ electrode + '/' + stim_type + '/' + amplitude + '/0' + electrode + stim_type + '_' + amplitude + '_' + network + '/spikes.csv')
                spikes_bg.append('exp' + exp + '/output/bkg/bkg_' + network + '/spikes.csv')
            elec_loc.append('../bio_components/stimulations/exp' + exp + '/elec_loc/' + electrode + '.csv')

    return dict(nodes_dirs=nodes, spikes_dirs=spikes, spikes_bg_dirs=spikes_bg, electrodes_dirs=elec_loc, radius=radius)

kwargs = get_params(1,3,'-',10,1)
print(kwargs)
spikes = get_spikes(**kwargs)
print(spikes)

def print_args(names=True, **kwargs):
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