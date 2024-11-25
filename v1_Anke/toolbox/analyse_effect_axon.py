import pandas as pd
from hdf5 import HDF5
import numpy as np
import os
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt

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
    z_in_between_min = 77
    z_in_between_max = 137

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

def projected_kernel_density_estimate(node_pos,n_spikes):
        '''
        Projection of the data before the Kernel Density Estimate is useful when trying the understand the spatial
        distribution of cell activity along a specific axis. The density estimate now reflects the distribution
        of the data along the projected data.
        If projection would only take place after the Kernel Density Estimate, then you just integrate the density
        function, but the density estimate is still based on a 2D-distribution.
        '''
        node_pos= node_pos[:,1:] #select only the y and z coordinates
        projected_points_z=node_pos[:,1]
        projected_points_y=node_pos[:,0]
        # Define grid along the projected axis
        grid_size = 100
        grid_z = np.linspace(min(projected_points_z), max(projected_points_z), grid_size).reshape(-1, 1)
        grid_y = np.linspace(min(projected_points_y), max(projected_points_y), grid_size).reshape(-1, 1)
        #print("projected points min and max z", min(projected_points_z), max(projected_points_z))
        #print("projected points min and max y", min(projected_points_y), max(projected_points_y))

        # Perform kernel density estimation
        kde = KernelDensity(bandwidth=200, kernel='gaussian') 

        kde_z=kde.fit(projected_points_z.reshape(-1, 1), sample_weight=n_spikes)
        density_z = np.exp(kde_z.score_samples(grid_z))
        kde_y=kde.fit(projected_points_y.reshape(-1, 1),sample_weight=n_spikes)
        density_y = np.exp(kde_y.score_samples(grid_y))

        fig, ax = plt.subplots(2, 1)
        # Plot 1D density function along the z axis 
        ax[0].plot(grid_z, density_z, color='red', linestyle='-')
        ax[0].set_xlabel('Distance along the z axis')
        ax[0].set_ylabel('Density')
        ax[0].set_title('1D Kernel Density Estimate along z axis')

        # Plot 1D density function along the z axis 
        ax[1].plot(grid_y, density_y, color='red', linestyle='-')
        ax[1].set_xlabel('Distance along the y axis')
        ax[1].set_ylabel('Density')
        ax[1].set_title('1D Kernel Density Estimate along y axis')
        plt.tight_layout()
        #plt.show()
        plt.close()

        return grid_y, grid_z, density_y, density_z

def plot1_kde(node_pos, n_spikes, pattern, mouse, amplitude):
    grid_y, grid_z, density_y, density_z = projected_kernel_density_estimate(node_pos, n_spikes)
    max_y_axis=grid_y[np.argmax(density_y)][0]
    max_z_axis=grid_z[np.argmax(density_z)][0]

    node_pos= node_pos[:,1:]
    max_spikes=np.max(n_spikes)
    n_spikes_norm=n_spikes/max_spikes

    electrode_0_zy=[16,300]
    electrode_1_zy=[198,300]
    electrode_2_zy=[16,170]
    electrode_3_zy=[198,170]

    fig = plt.figure(figsize=(8,12))

    if pattern==0:
        plt.scatter(electrode_1_zy[0], electrode_1_zy[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)
    elif pattern==4:
        plt.scatter(electrode_1_zy[0], electrode_1_zy[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)
    elif pattern==5:
        plt.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
    elif pattern==6:
        plt.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
    elif pattern==8:
        plt.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
        plt.scatter(electrode_3_zy[0], electrode_3_zy[1], color='yellow', s=110, marker='s', label='Return electrode 2 in L2/3', zorder=3)
    else:
        plt.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
        plt.scatter(electrode_3_zy[0], electrode_3_zy[1], color='yellow', s=110, marker='s', label='Return electrode 2 in L2/3', zorder=3)

    #plt.axline(electrode_0_zy, electrode_1_zy, color='limegreen', label='Along layer')
    #plt.axline(electrode_0_zy, electrode_2_zy, color='darkgreen', label='Along column')
    plt.scatter(node_pos[:,1], node_pos[:,0], s=90, c="blue", alpha=n_spikes_norm)
    plt.scatter(electrode_0_zy[0], electrode_0_zy[1], color='orange', s=110, marker='s', label='Central electrode', zorder=3)
    #plt.scatter(max_z_axis, electrode_0_zy[1], color='red', marker='*', s=120, label='Max density', zorder=3)
    #plt.scatter(electrode_0_zy[0], max_y_axis, color='red', marker='*', s=120, zorder=3)
    #plt.scatter(max_z_axis,max_y_axis, color='red', marker='*', s=120, zorder=3)

    y_min = 100 # border layer 1 and 2/3
    y_max = 430 # border layer 4 and 5
    z_electrode_min = -14
    z_electrode_max = 46
    z_in_between_min = 77
    z_in_between_max = 137

    plt.plot([z_electrode_min, z_electrode_max], [y_min, y_min], color='r')  # Top horizontal line
    plt.plot([z_electrode_min, z_electrode_max], [y_max, y_max], color='r')  # Bottom horizontal line
    plt.plot([z_electrode_min, z_electrode_min], [y_min, y_max], color='r')  # Left vertical line
    plt.plot([z_electrode_max, z_electrode_max], [y_min, y_max], color='r')  # Right vertical line

    plt.plot([z_in_between_min, z_in_between_max], [y_min, y_min], color='r')  # Top horizontal line
    plt.plot([z_in_between_min, z_in_between_max], [y_max, y_max], color='r')  # Bottom horizontal line
    plt.plot([z_in_between_min, z_in_between_min], [y_min, y_max], color='r')  # Left vertical line
    plt.plot([z_in_between_max, z_in_between_max], [y_min, y_max], color='r')  # Right vertical line



    plt.xlabel('Z Coordinate')
    plt.ylabel('Y Coordinate')
    plt.xlim([-400, 400])
    plt.ylim([0, 800])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    plt.legend(fontsize='12', loc='upper right')

    pattern_title="Parallel to cortical columns. Pattern"+str(pattern)+". M"+str(mouse)+". Amplitude "+ str(amplitude)+"."
    plt.title(pattern_title)
    #plt.savefig('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_4/plots_column/exp_4_density_1d_kde_p'+str(pattern)+'_m_'+str(mouse)+'a_'+str(amplitude)+'.png')
    #plt.close()
    plt.show()
    return max_y_axis, max_z_axis

def plot1_kde_double(node_pos_A, n_spikes_A, node_pos_B, n_spikes_B, pattern_A, mouse_A, amplitude_A, pattern_B, mouse_B, amplitude_B):
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))

    # Plot for A
    grid_y_A, grid_z_A, density_y_A, density_z_A = projected_kernel_density_estimate(node_pos_A, n_spikes_A)
    max_y_axis_A = grid_y_A[np.argmax(density_y_A)][0]
    max_z_axis_A = grid_z_A[np.argmax(density_z_A)][0]

    node_pos_A = node_pos_A[:, 1:]
    max_spikes_A = np.max(n_spikes_A)
    n_spikes_norm_A = n_spikes_A / max_spikes_A

    electrode_0_zy_A = [16, 300]
    electrode_1_zy_A = [198, 300]
    electrode_2_zy_A = [16, 170]
    electrode_3_zy_A = [198, 170]

    axs[0].scatter(electrode_1_zy_A[0], electrode_1_zy_A[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)

    # Plot positions and spikes of A
    axs[0].scatter(node_pos_A[:, 1], node_pos_A[:, 0], s=90, c="blue", alpha=n_spikes_norm_A)
    axs[0].scatter(electrode_0_zy_A[0], electrode_0_zy_A[1], color='orange', s=110, marker='s', label='Central electrode', zorder=3)
    #axs[0].scatter(max_z_axis_A, electrode_0_zy_A[1], color='red', marker='*', s=120, label='Max density', zorder=3)
    #axs[0].scatter(electrode_0_zy_A[0], max_y_axis_A, color='red', marker='*', s=120, zorder=3)
    #axs[0].scatter(max_z_axis_A, max_y_axis_A, color='red', marker='*', s=120, zorder=3)

    y_min_A = 100  # border layer 1 and 2/3
    y_max_A = 430  # border layer 4 and 5
    z_electrode_min_A = -14
    z_electrode_max_A = 46
    z_in_between_min_A = 77
    z_in_between_max_A = 137

    # Plot box for A
    axs[0].plot([z_electrode_min_A, z_electrode_max_A], [y_min_A, y_min_A], color='r')  # Top horizontal line
    axs[0].plot([z_electrode_min_A, z_electrode_max_A], [y_max_A, y_max_A], color='r')  # Bottom horizontal line
    axs[0].plot([z_electrode_min_A, z_electrode_min_A], [y_min_A, y_max_A], color='r')  # Left vertical line
    axs[0].plot([z_electrode_max_A, z_electrode_max_A], [y_min_A, y_max_A], color='r')  # Right vertical line
    axs[0].plot([z_in_between_min_A, z_in_between_max_A], [y_min_A, y_min_A], color='r')  # Top horizontal line
    axs[0].plot([z_in_between_min_A, z_in_between_max_A], [y_max_A, y_max_A], color='r')  # Bottom horizontal line
    axs[0].plot([z_in_between_min_A, z_in_between_min_A], [y_min_A, y_max_A], color='r')  # Left vertical line
    axs[0].plot([z_in_between_max_A, z_in_between_max_A], [y_min_A, y_max_A], color='r')  # Right vertical line

    axs[0].set_xlabel('Z Coordinate')
    axs[0].set_ylabel('Y Coordinate')
    axs[0].set_xlim([-400, 400])
    axs[0].set_ylim([0, 800])
    axs[0].invert_xaxis()
    axs[0].invert_yaxis()
    axs[0].set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    axs[0].legend(fontsize='12', loc='upper right')

    pattern_title_A = "Parallel to cortical columns. Original. Pattern" + str(pattern_A) + ". M" + str(mouse_A) + ". Amplitude " + str(amplitude_A) + "."
    axs[0].set_title(pattern_title_A)

    # Plot for B
    grid_y_B, grid_z_B, density_y_B, density_z_B = projected_kernel_density_estimate(node_pos_B, n_spikes_B)
    max_y_axis_B = grid_y_B[np.argmax(density_y_B)][0]
    max_z_axis_B = grid_z_B[np.argmax(density_z_B)][0]

    node_pos_B = node_pos_B[:, 1:]
    max_spikes_B = np.max(n_spikes_B)
    n_spikes_norm_B = n_spikes_B / max_spikes_B

    electrode_0_zy_B = [16, 300]
    electrode_1_zy_B = [198, 300]
    electrode_2_zy_B = [16, 170]
    electrode_3_zy_B = [198, 170]

    axs[1].scatter(electrode_1_zy_B[0], electrode_1_zy_B[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)

    # Plot positions and spikes of B
    axs[1].scatter(node_pos_B[:, 1], node_pos_B[:, 0], s=90, c="blue", alpha=n_spikes_norm_B)
    axs[1].scatter(electrode_0_zy_B[0], electrode_0_zy_B[1], color='orange', s=110, marker='s', label='Central electrode', zorder=3)
    #axs[1].scatter(max_z_axis_B, electrode_0_zy_B[1], color='red', marker='*', s=120, label='Max density', zorder=3)
    #axs[1].scatter(electrode_0_zy_B[0], max_y_axis_B, color='red', marker='*', s=120, zorder=3)
    #axs[1].scatter(max_z_axis_B, max_y_axis_B, color='red', marker='*', s=120, zorder=3)

    axs[1].plot([z_electrode_min_A, z_electrode_max_A], [y_min_A, y_min_A], color='r')  # Top horizontal line
    axs[1].plot([z_electrode_min_A, z_electrode_max_A], [y_max_A, y_max_A], color='r')  # Bottom horizontal line
    axs[1].plot([z_electrode_min_A, z_electrode_min_A], [y_min_A, y_max_A], color='r')  # Left vertical line
    axs[1].plot([z_electrode_max_A, z_electrode_max_A], [y_min_A, y_max_A], color='r')  # Right vertical line
    axs[1].plot([z_in_between_min_A, z_in_between_max_A], [y_min_A, y_min_A], color='r')  # Top horizontal line
    axs[1].plot([z_in_between_min_A, z_in_between_max_A], [y_max_A, y_max_A], color='r')  # Bottom horizontal line
    axs[1].plot([z_in_between_min_A, z_in_between_min_A], [y_min_A, y_max_A], color='r')  # Left vertical line
    axs[1].plot([z_in_between_max_A, z_in_between_max_A], [y_min_A, y_max_A], color='r')  # Right vertical line

    axs[1].set_xlabel('Z Coordinate')
    axs[1].set_ylabel('Y Coordinate')
    axs[1].set_xlim([-400, 400])
    axs[1].set_ylim([0, 800])
    axs[1].invert_xaxis()
    axs[1].invert_yaxis()
    axs[1].set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    axs[1].legend(fontsize='12', loc='upper right')

    pattern_title_B = "Parallel to cortical columns. Active axons. Pattern " + str(pattern_A) + ". M" + str(mouse_A) + ". Amplitude " + str(amplitude_A) + "."
    axs[1].set_title(pattern_title_B)

    #plt.savefig('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_4/compare_exp_2/comparison_exp2_and_4.png')
    plt.show()
    return

exp_A=2
pattern_A=0
mouse_A=0
amplitude_A=10
node_pos_A, n_spikes_A = get_spikes(exp=exp_A,pattern=pattern_A,mouse=mouse_A,amplitude=amplitude_A)

exp_B=4
pattern_B=0
mouse_B=0
amplitude_B=10
node_pos_B, n_spikes_B = get_spikes(exp=exp_B,pattern=pattern_B,mouse=mouse_B,amplitude=amplitude_B)

positions_filtered_A, spikes_filtered_A, threshold_A = filter_spikes(node_pos_A, n_spikes_A)
positions_filtered_B, spikes_filtered_B, threshold_B = filter_spikes(node_pos_B, n_spikes_B)

#in_between_norm_A= densities_active_neurons(positions_filtered_A, spikes_filtered_A)
#n_between_norm_B= densities_active_neurons(positions_filtered_B, spikes_filtered_B)

#max_y_axis_A, max_z_axis_A = plot1_kde(positions_filtered_A, spikes_filtered_A, pattern_A, mouse_A,amplitude_A)
plot1_kde_double(positions_filtered_A, spikes_filtered_A, positions_filtered_B, spikes_filtered_B, pattern_A, mouse_A, amplitude_A, pattern_B, mouse_B, amplitude_B)