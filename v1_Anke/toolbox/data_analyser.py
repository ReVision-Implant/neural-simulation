import pandas as pd
from hdf5 import HDF5
import numpy as np
import os

def get_spikes(nodes_dirs, spikes_dirs, spikes_bkg_dirs, radius=None, depth=None, v1=True, **kwargs):
    """Get spikes and node positions from network and output files.

    :param nodes_dirs: directories that point to network/nodes.h5
    :type nodes_dirs: path or list thereof
    :param spikes_dirs: directories that point to output/spikes.csv
    :type spikes_dirs: path or list thereof
    :param spikes_bkg_dirs: directories that point to output/bkg/spikes.csv
    :type spikes_bkg_dirs: path or list thereof
    :param v1: defaults to True.
    :type v1: bool, optional
    :return: node positions
    :rtype: ndarray
    """    
        
    nodes_dirs = [nodes_dirs] if not isinstance(nodes_dirs, list) else nodes_dirs
    spikes_dirs = [spikes_dirs] if not isinstance(spikes_dirs, list) else spikes_dirs
    spikes_bkg_dirs = [spikes_bkg_dirs] if not isinstance(spikes_bkg_dirs, list) else spikes_bkg_dirs

    assert len(nodes_dirs) == len(spikes_dirs) == len(spikes_bkg_dirs)

    node_pos = np.zeros((1,3))
    n_spikes = np.zeros((1,1)) 

    for i in range(len(nodes_dirs)):

        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        spikes_bkg_dir = spikes_bkg_dirs[i]

        node_pos_temp = HDF5(nodes_dir, v1=v1).positions

        n_spikes_temp = np.zeros(np.shape(node_pos_temp)[0])

        spikes = pd.read_csv(spikes_dir, sep='\s+')
        for ind in spikes.index:
            if spikes['timestamps'][ind] < 100:
                n_spikes_temp[spikes['node_ids'][ind]] += 1

        if spikes_bkg_dirs is not None:
            spikes_bkg_dir = spikes_bkg_dirs[0] if isinstance(spikes_bkg_dirs, list) else spikes_bkg_dirs 
            spikes_bkg = pd.read_csv(spikes_bkg_dir, sep='\s+')
            for ind in spikes_bkg.index:
                if spikes_bkg['timestamps'][ind] < 100:
                    n_spikes_temp[spikes_bkg['node_ids'][ind]] = max(0, n_spikes_temp[spikes_bkg['node_ids'][ind]] - 1)

        node_pos = np.vstack((node_pos, node_pos_temp))
        n_spikes = np.append(n_spikes, n_spikes_temp)

    return node_pos, n_spikes

path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
node_dirs_A = [path+'/virtual_mice_mask/mouse_0/v1_nodes.h5']
spike_dirs_A = [path+'/exp_2/output/pattern_0/amplitude_10/mouse_0/spikes.csv']
spike_bkg_dirs_A= [path+'/exp_2/output/bkg/mouse_0/spikes.csv']
node_pos_A, n_spikes_A = get_spikes(nodes_dirs = node_dirs_A, spikes_dirs = spike_dirs_A, spikes_bkg_dirs = spike_bkg_dirs_A)

node_dirs_B = [path+'/virtual_mice_mask/mouse_0/v1_nodes.h5']
spike_dirs_B = [path+'/exp_2/output/pattern_4/amplitude_10/mouse_0/spikes.csv']
spike_bkg_dirs_B= [path+'/exp_2/output/bkg/mouse_0/spikes.csv']
node_pos_B, n_spikes_B = get_spikes(nodes_dirs = node_dirs_B, spikes_dirs = spike_dirs_B, spikes_bkg_dirs = spike_bkg_dirs_B)

#print(n_spikes_A.shape)
#print(n_spikes_B.shape)
#print(np.max(n_spikes_A))
#print(np.max(n_spikes_B))
