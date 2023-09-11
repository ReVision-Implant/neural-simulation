import pandas as pd
from hdf5 import HDF5
import numpy as np
import os

def get_spikes(nodes_dirs, spikes_dirs, spikes_bkg_dirs, radius=None, depth=None, v1=True, **kwargs):
    
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
        
        if radius is not None:
            circle = node_pos_temp[:,0]**2 + node_pos_temp[:,2]**2 < radius**2
            node_pos_temp = node_pos_temp[circle,:]
            n_spikes_temp = n_spikes_temp[circle]

        if depth is not None:
            layer = (node_pos_temp[:,1] >= depth[0])*(node_pos_temp[:,1] < depth[1])  
            node_pos_temp = node_pos_temp[layer,:]
            n_spikes_temp = n_spikes_temp[layer]

        node_pos = np.vstack((node_pos, node_pos_temp))
        n_spikes = np.append(n_spikes, n_spikes_temp)

    return node_pos, n_spikes

def get_grid(node_pos, n_spikes, radius, **kwargs):

    grid = np.zeros((radius*2,radius*2))
    for node in range(node_pos.shape[0]):
        grid_el = (np.floor(node_pos[node,[0,1,2]] + radius)).astype(int)
        grid[grid_el[2],grid_el[0]] += n_spikes[node]

    return grid

def get_centroid_cov(node_pos, n_spikes, **kwargs):

    # Get centroid and cov matrix from spikes+locations from previous step
    centroid = np.average(node_pos, axis=0, weights=n_spikes)
    cov = np.cov(node_pos.T, fweights=n_spikes)
    
    return centroid, cov


if __name__ == '__file__':
    None