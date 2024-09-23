import pandas as pd
from hdf5 import HDF5
import numpy as np
import os
#from scipy.stats import wilcoxon
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import math


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

    path ='C:/Users/ankev/OneDrive/Documenten/Github/ReVision/neural-simulation/v1_Anke'
    
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
    #print("average", avg_spikes)
    std_spikes = np.std(n_spikes)
    #print("standard dev", std_spikes)    

    threshold = avg_spikes + 3*std_spikes
    #print("threshold", threshold)
    n_spikes_filtered=[]
    filtered_indices=[]
    for index, value in enumerate(n_spikes):
        if value >= threshold or value >=threshold:
            n_spikes_filtered.append(value)
            filtered_indices.append(index)

    #print("before filtering node pos shape:", node_pos.shape)            
    node_pos_filtered = node_pos[filtered_indices]
    #print("after filtering node pos shape:", node_pos_filtered.shape)   
    n_spikes_filtered= np.array(n_spikes_filtered)
    #print("after filtering n_ spikes shape:", n_spikes_filtered.shape) 

    return node_pos_filtered, n_spikes_filtered, threshold

def correlation(n_spikes_A, n_spikes_B,pattern_A,pattern_B,threshold_A, threshold_B):
        '''
        Use the Spearman correlation test to get a p-value as index of separability between the two neuronal populations.
        '''
        n_spikes_A_filtered=[]
        n_spikes_B_filtered=[]
        for value1, value2 in zip(n_spikes_A, n_spikes_B):
            if value1 >= threshold_A or value2 >=threshold_B:
            #print(value1,value2)
                n_spikes_A_filtered.append(value1)
                n_spikes_B_filtered.append(value2)

        n_spikes_A= n_spikes_A_filtered
        n_spikes_B= n_spikes_B_filtered

        statistic, pvalue = spearmanr(n_spikes_A, n_spikes_B, alternative ="greater")
        
        # Calculate the confidence interval
        n = len(n_spikes_A)
        stderr = 1.0 / math.sqrt(n - 3)
        z_score = 1.96  # For a 95% confidence interval
        delta = z_score * stderr
        lower_bound = math.tanh(math.atanh(statistic) - delta)
        upper_bound = math.tanh(math.atanh(statistic) + delta)
    
        print('The Spearman correlation coefficient for stim patterns ' + str(pattern_A) + ' and ' + str(pattern_B) +
          ' is ' + str(round(statistic, 2)) + ', the p-value is ' + str(round(pvalue, 4)) +
          ', and the 95% confidence interval is [' + str(round(lower_bound, 2)) + ', ' + str(round(upper_bound, 2)) + ']')

        return(statistic, pvalue)

def overlap(n_spikes_A, n_spikes_B, threshold_A, threshold_B): 
        
        n_spikes_A_filtered=[]
        n_spikes_B_filtered=[]
        activity_A =[]
        activity_B =[]
        for value1, value2 in zip(n_spikes_A, n_spikes_B):
            if value1 >= threshold_A or value2 >=threshold_B:
            #print(value1,value2)
                n_spikes_A_filtered.append(value1)
                n_spikes_B_filtered.append(value2)

                if value1 >= threshold_A:
                     activity_A.append(1)
                else:
                     activity_A.append(0)
                if value1 >= threshold_B:
                     activity_B.append(1)
                else:
                     activity_B.append(0)     

        n_spikes_A= n_spikes_A_filtered
        n_spikes_B= n_spikes_B_filtered

        #print(threshold_A)
        #print(n_spikes_A)
        #print(activity_A)

        overlap = [1 if (activity_A[i] == 1 and activity_B[i] == 1) else 0 for i in range(len(activity_A))]

        active_neurons_total = len(activity_A)
        active_neurons_A = activity_A.count(1)
        active_neurons_B = activity_B.count(1)
        active_neurons_overlap = overlap.count(1)

        if active_neurons_total > 0:
            return np.round(active_neurons_overlap/active_neurons_total,2), active_neurons_total
        else:
            return 0, 0

def get_electrode_angles(central, elec1, elec2):
        P12 = np.sqrt((central[2]-elec1[2])**2 + (central[1]-elec1[1])**2)
        P13 = np.sqrt((central[2]-elec2[2])**2 + (central[1]-elec2[1])**2)
        P23 = np.sqrt((elec1[2]-elec2[2])**2 + (elec1[1]-elec2[1])**2)
        return np.rad2deg(np.arccos((P12**2 + P13**2 - P23**2)/(2*P12*P13)))

def correlation_per_angle(exp=[4,5], patterns=[0,1], mouse=[0,0], amplitude=[10,10]): #hier aan gewerkt
        angles = []
        correlations = []
        overlaps = []
    
        for i in range(len(patterns)):
            node_pos_i, n_spikes_i= get_spikes(exp[i],patterns[i],mouse[i],amplitude[i])
            pos_filtered_i, spikes_filtered_i, threshold_i = filter_spikes(node_pos_i,n_spikes_i)
            
            for j in range(i+1,len(patterns)):
                node_pos_j, n_spikes_j= get_spikes(exp[j],patterns[j],mouse[j],amplitude[j])
                pos_filtered_j, spikes_filtered_j, threshold_j = filter_spikes(node_pos_j,n_spikes_j)

                #correlations
                statistic, pvalue = correlation(n_spikes_i, n_spikes_j,patterns[i],patterns[j],threshold_i, threshold_j)
                correlations.append(statistic)

                #overlaps
                overlaps.append(overlap(n_spikes_i, n_spikes_j, threshold_i, threshold_j))

                #angles
                el_coordinates = {
                    0: [-9, 300, 16],
                    1: [-9, 300, 198],
                    2: [-9, 170, 16],
                    3: [-9, 170, 380],
                    4: [-9, 170, 380],
                    5: [-9, 300, 380],
                    6: [-9, 170, -166],
                    7: [-9, 170, -348],
                    8: [-9, 300, -348]
                }

                if exp[i]==4:
                    exp_i = {
                        0: [0,1],
                        5: [0,2],
                        #7: [0,2,3]
                        }
                else:
                    exp_i = {
                        0: [0,5],
                        1: [0,4],
                        #2: [0,3,4],
                        #3: [0,3,4],
                        4: [0,7],
                        5: [0,8],
                        6: [0,9],
                        7: [0,2,6],
                        #8: [0,6,7],
                        9: [0,3]
                        } 

                if exp[j]==4:
                    exp_j = {
                        0: [0,1],
                        5: [0,2],
                        #7: [0,2,3]
                        }
                else:
                    exp_j = {
                        0: [0,5],
                        1: [0,4],
                        #2: [0,3,4],
                        #3: [0,3,4],
                        4: [0,7],
                        5: [0,8],
                        6: [0,9],
                        7: [0,2,6],
                        #8: [0,6,7],
                        9: [0,3]
                        } 
                       
            
                electrodes_i = exp_i[patterns[i]]
                electrodes_j = exp_i[patterns[j]]
                return_el_i = electrodes_i[1]
                #print("return electrode i = ", return_el_i)
                return_el_j = electrodes_j[1]
                #print("return electrode j = ", return_el_j)

                location_central_el = [-9, 300, 16]
                location_return_i = el_coordinates[return_el_i]
                #print("location electrode i = ", location_return_i)
                location_return_j = el_coordinates[return_el_j]
                #print("location electrode j = ", location_return_j)


                angles.append(get_electrode_angles(location_central_el, location_return_i,location_return_j))
       
        angles, correlations, overlaps = zip(*sorted(zip(angles, correlations, overlaps))) # Sort the correlations in ascending order
        return angles, correlations, overlaps

##angles klaar, figuur proberen
# daarna symmetry plot
# centroid and depth plots

##test the code
exp_A=4
pattern_A=0
mouse_A=0
amplitude_A=10
#node_pos_A, n_spikes_A = get_spikes(exp=exp_A,pattern=pattern_A,mouse=mouse_A,amplitude=amplitude_A)

exp_B = 4
pattern_B=5
mouse_B=0
amplitude_B=10
#node_pos_B, n_spikes_B = get_spikes(exp=exp_B,pattern=pattern_B,mouse=mouse_B,amplitude=amplitude_B)

#positions_filtered_A, spikes_filtered_A, threshold_A = filter_spikes(node_pos_A, n_spikes_A)
#positions_filtered_B, spikes_filtered_B, threshold_B = filter_spikes(node_pos_B, n_spikes_B)

#overlap(n_spikes_A, n_spikes_B, threshold_A, threshold_B)

correlation_per_angle(exp=[4,4], patterns=[0,5], mouse=[0,0], amplitude=[10,10])