import pandas as pd
from hdf5 import HDF5
import numpy as np
import os


import numpy as np
import pandas as pd
from h5py import File as HDF5

def get_spikes(nodes_dirs, spikes_dirs, spikes_bkg_dirs, **kwargs):
    """Get spikes and node positions from network and output files.

    :param nodes_dirs: directories that point to network/nodes.h5
    :type nodes_dirs: path or list thereof
    :param spikes_dirs: directories that point to output/spikes.csv
    :type spikes_dirs: path or list thereof
    :param spikes_bkg_dirs: directories that point to output/bkg/spikes.csv
    :type spikes_bkg_dirs: path or list thereof
    :param v1: defaults to True.
    :type v1: bool, optional
    :return: average node positions and average number of spikes for each node
    :rtype: ndarray
    """    
        
    nodes_dirs = [nodes_dirs] if not isinstance(nodes_dirs, list) else nodes_dirs
    spikes_dirs = [spikes_dirs] if not isinstance(spikes_dirs, list) else spikes_dirs
    spikes_bkg_dirs = [spikes_bkg_dirs] if not isinstance(spikes_bkg_dirs, list) else spikes_bkg_dirs

    assert len(nodes_dirs) == len(spikes_dirs) == len(spikes_bkg_dirs)

    node_pos_total = np.zeros((0,3)) # Initialize total node positions
    n_spikes_total = np.zeros((0,1)) # Initialize total spikes count

    for i in range(len(nodes_dirs)):
        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        spikes_bkg_dir = spikes_bkg_dirs[i]

        node_pos_temp = HDF5(nodes_dir, v1=True).positions
        n_spikes_temp = np.zeros(np.shape(node_pos_temp)[0])

        spikes = pd.read_csv(spikes_dir, sep='\s+')
        for ind in spikes.index:
            if spikes['timestamps'][ind] < 100:
                n_spikes_temp[spikes['node_ids'][ind]] += 1

        if spikes_bkg_dirs is not None:
            spikes_bkg = pd.read_csv(spikes_bkg_dir, sep='\s+')
            for ind in spikes_bkg.index:
                if spikes_bkg['timestamps'][ind] < 100:
                    n_spikes_temp[spikes_bkg['node_ids'][ind]] = max(0, n_spikes_temp[spikes_bkg['node_ids'][ind]] - 1)

        node_pos_total = np.vstack((node_pos_total, node_pos_temp))
        n_spikes_total = np.append(n_spikes_total, n_spikes_temp)

    # Calculate averages
    avg_node_pos = node_pos_total.mean(axis=0)
    avg_n_spikes = n_spikes_total.mean(axis=0)

    return avg_node_pos, avg_n_spikes


path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
node_dirs_0 = [path+'/virtual_mice_mask/mouse_0/v1_nodes.h5',path+'/virtual_mice_mask/mouse_1/v1_nodes.h5']
spike_dirs_0 = [path+'/exp_2/output/pattern_0/amplitude_10/mouse_0/spikes.csv',path+'/exp_2/output/pattern_0/amplitude_10/mouse_1/spikes.csv']
spike_bkg_dirs_0= [path+'/exp_2/output/bkg/mouse_0/spikes.csv', path+'/exp_2/output/bkg/mouse_1/spikes.csv']
node_positions_pattern_0_amp_10, n_spikes_pattern_0_amp_10 = get_spikes(nodes_dirs = node_dirs_0, spikes_dirs = spike_dirs_0, spikes_bkg_dirs = spike_bkg_dirs_0)

print(node_positions_pattern_0_amp_10[147996])
print(n_spikes_pattern_0_amp_10[147996])
max_value = np.max(n_spikes_pattern_0_amp_10)
#print(max_value)