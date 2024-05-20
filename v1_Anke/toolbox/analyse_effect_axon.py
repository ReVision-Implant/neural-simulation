import pandas as pd
from hdf5 import HDF5
import numpy as np
import os

def get_spikes(exp,pattern,mouse,amplitude, v1=True, **kwargs):
    """Get spikes and node positions from network and output files.

    :param nodes_dirs: directories that point to network/nodes.h5
    :type nodes_dirs: path or list thereof
    :param spikes_dirs: directories that point to output/spikes.csv
    :type spikes_dirs: path or list thereof
    :param v1: defaults to True.
    :type v1: bool, optional
    :return: node positions
    :rtype: ndarray
    """    

    path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
    
    nodes_dirs= [str(path)+'/virtual_mice_mask/mouse_'+str(mouse)+'/v1_nodes.h5']
    spikes_dirs= [str(path)+'/exp_'+str(exp)+'/output/pattern_'+str(pattern)+'/amplitude_'+str(amplitude)+'/mouse_'+str(mouse)+'/spikes.csv']
    #spikes_bkg_dirs= [str(path)+'/exp_'+str(exp)+'/output/bkg/mouse_'+str(mouse)+'/spikes.csv']
        
    nodes_dirs = [nodes_dirs] if not isinstance(nodes_dirs, list) else nodes_dirs
    spikes_dirs = [spikes_dirs] if not isinstance(spikes_dirs, list) else spikes_dirs
    #spikes_bkg_dirs = [spikes_bkg_dirs] if not isinstance(spikes_bkg_dirs, list) else spikes_bkg_dirs

    #assert len(nodes_dirs) == len(spikes_dirs) == len(spikes_bkg_dirs)
    assert len(nodes_dirs) == len(spikes_dirs)

    node_pos = np.zeros((1,3))
    n_spikes = np.zeros((1,1)) 

    for i in range(len(nodes_dirs)):

        nodes_dir = nodes_dirs[i]
        spikes_dir = spikes_dirs[i]
        #spikes_bkg_dir = spikes_bkg_dirs[i]

        node_pos_temp = HDF5(nodes_dir, v1=v1).positions

        n_spikes_temp = np.zeros(np.shape(node_pos_temp)[0])

        spikes = pd.read_csv(spikes_dir, sep='\s+')
        for ind in spikes.index:
            if spikes['timestamps'][ind] < 100:
                n_spikes_temp[spikes['node_ids'][ind]] += 1

        #if spikes_bkg_dirs is not None:
        #    spikes_bkg_dir = spikes_bkg_dirs[0] if isinstance(spikes_bkg_dirs, list) else spikes_bkg_dirs 
        #    spikes_bkg = pd.read_csv(spikes_bkg_dir, sep='\s+')
        #    for ind in spikes_bkg.index:
        #        if spikes_bkg['timestamps'][ind] < 100:
         #           n_spikes_temp[spikes_bkg['node_ids'][ind]] = max(0, n_spikes_temp[spikes_bkg['node_ids'][ind]] - 1)

        node_pos = np.vstack((node_pos, node_pos_temp))
        n_spikes = np.append(n_spikes, n_spikes_temp)

    return node_pos, n_spikes

def filter_spikes(node_pos, n_spikes):
    non_zero_indices = np.nonzero(n_spikes)
    node_pos = node_pos[non_zero_indices]
    n_spikes= n_spikes[non_zero_indices]

    avg_spikes = np.mean(n_spikes)
    print("average", avg_spikes)
    std_spikes = np.std(n_spikes)
    print("standard dev", std_spikes)    

    threshold = avg_spikes + 3*std_spikes
    print("threshold", threshold)
    n_spikes_filtered=[]
    filtered_indices=[]
    for index, value in enumerate(n_spikes):
        if value >= threshold or value >=threshold:
            n_spikes_filtered.append(value)
            filtered_indices.append(index)

    print("before filtering node pos shape:", node_pos.shape)            
    node_pos_filtered = node_pos[filtered_indices]
    print("after filtering node pos shape:", node_pos_filtered.shape)   
    n_spikes_filtered= np.array(n_spikes_filtered)
    print("after filtering n_ spikes shape:", n_spikes_filtered.shape) 

    return node_pos_filtered, n_spikes_filtered, threshold


def densities_active_neurons(node_pos_filtered, n_spikes_filtered):
    # Extracting y and z coordinates from node_pos_filtered
    node_pos = node_pos_filtered[:, 1:]

    # Define the ranges
    y_min = 100 # border layer 1 and 2/3
    y_max = 430 # border layer 4 and 5
    z_electrode_min = -14
    z_electrode_max = 46
    z_in_between_min = 61
    z_in_between_max = 121

    # Initialize variables to store summed spike rates
    around_electrode = 0
    in_between = 0
    number_neurons_around_electrode = 0
    number_neurons_in_between = 0

    # Iterate through each node
    for i in range(len(node_pos)):
        y = node_pos[i, 0]
        z = node_pos[i, 1]

        # Check if the node falls within the ranges
        if (y_min <= y <= y_max) and (z_electrode_min <= z <= z_electrode_max):
            around_electrode += n_spikes_filtered[i]
            number_neurons_around_electrode+=1

        if (y_min <= y <= y_max) and (z_in_between_min <= z <= z_in_between_max):
            in_between += n_spikes_filtered[i]
            number_neurons_in_between+=1
    
    #normalize with the around electrode spike rate
    in_between_norm=in_between/around_electrode
    print("spikes around the electrode ", around_electrode)
    print("number of active neurons around electrode ", number_neurons_around_electrode)
    print("spikes in between the electrodes ", in_between)
    print("number of active neurons in between ", number_neurons_in_between)
    print("spikes in between the electrode divided by the spikes around the electrode ", in_between_norm)

    return in_between_norm

exp_A=4
pattern_A=0
mouse_A=0
amplitude_A=10
node_pos_A, n_spikes_A = get_spikes(exp=exp_A,pattern=pattern_A,mouse=mouse_A,amplitude=amplitude_A)

exp_B=2
pattern_B=0
mouse_B=0
amplitude_B=10
node_pos_B, n_spikes_B = get_spikes(exp=exp_B,pattern=pattern_B,mouse=mouse_B,amplitude=amplitude_B)

positions_filtered_A, spikes_filtered_A, threshold_A = filter_spikes(node_pos_A, n_spikes_A)
positions_filtered_B, spikes_filtered_B, threshold_B = filter_spikes(node_pos_B, n_spikes_B)

in_between_norm_A= densities_active_neurons(positions_filtered_A, spikes_filtered_A)
in_between_norm_B= densities_active_neurons(positions_filtered_B, spikes_filtered_B)