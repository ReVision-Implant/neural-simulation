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

#same as the other data_analyser file but now preprocessing with standard deviations instead of threshold

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

def Pearsoncorrel(n_spikes_A, n_spikes_B,pattern_A,pattern_B,threshold_A, threshold_B):
        '''
        Use the Pearson correlation test to get a p-value as index of separability between the two neuronal populations.
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

        statistic, pvalue = pearsonr(n_spikes_A, n_spikes_B, alternative ="greater")
        print('The Pearson correlation coefficient for stim patterns '+str(pattern_A)+' and '+str(pattern_B)+' is ' + str(round(statistic,2))+', the pvalue is '+str(round(pvalue,4)))
        
        plt.figure()
        plt.scatter(n_spikes_A, n_spikes_B, s=30)
        plt.xlabel('Spike rates for pattern ' + str(pattern_A))
        plt.ylabel('Spike rates for pattern ' + str(pattern_B))
        plt.title("background substracted + threshold activity = average + 3SD")
        plt.show()
    

        return(statistic, pvalue)

def Spearmancorr(n_spikes_A, n_spikes_B,pattern_A,pattern_B,threshold_A, threshold_B):
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
        
        plt.figure()
        plt.scatter(n_spikes_A, n_spikes_B, s=30)
        plt.xlabel('Spike rates for pattern ' + str(pattern_A))
        plt.ylabel('Spike rates for pattern ' + str(pattern_B))
        plt.title("background substracted + threshold activity = average + 3SD")
        #plt.show()
    

        return(statistic, pvalue)

def kernel_density_estimate(node_pos, n_spikes, pattern):
        '''
        2D Kernel Density Estimate of the data
        '''
        node_pos= node_pos[:,1:] #select only the y and z coordinates
        kde = KernelDensity(bandwidth=200, kernel='gaussian') # Choose model and parameters
        ###vanaf hier verder werken
        kde.fit(node_pos, sample_weight=n_spikes) # Train model

        grid_size = 100 # 100 points in x and in y direction
        y_grid, z_grid = np.meshgrid(np.linspace(0, 800, grid_size), np.linspace(-400, 400, grid_size))
        grid_points = np.vstack([y_grid.ravel(), z_grid.ravel()]).T
        density = np.exp(kde.score_samples(grid_points)).reshape(y_grid.shape) # Evaluate model for all points on the grid

        fig = plt.figure()
        plt.pcolormesh(z_grid, y_grid, density, shading='auto')
        plt.scatter(node_pos[:,1], node_pos[:,0], c=n_spikes, cmap='viridis', edgecolors='k', linewidths=1)
        plt.colorbar(label='Values')
        plt.xlabel('Z Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([-400, 400])
        plt.ylim([0, 800])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        # plt.legend()
        plt.title('Kernel Density Estimate for stim. pattern ' + str(pattern))
        #plt.show()
        plt.close()

        return y_grid, z_grid, density

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

def full_kde(node_pos, n_spikes, pattern, mouse, amplitude):
    grid_y_2D,grid_z_2D, density_2D=kernel_density_estimate(node_pos, n_spikes, pattern)
    grid_y, grid_z, density_y, density_z = projected_kernel_density_estimate(node_pos, n_spikes)
    max_y_axis=grid_y[np.argmax(density_y)][0]
    max_z_axis=grid_z[np.argmax(density_z)][0]

    node_pos= node_pos[:,1:]
    max_spikes=np.max(n_spikes)
    #print("max number spikes", max_spikes)
    n_spikes_norm=n_spikes/max_spikes
    #print(n_spikes_norm)

    electrode_0_zy=[16,300]
    electrode_1_zy=[198,300]
    electrode_2_zy=[16,170]
    electrode_3_zy=[198,170]

    fig = plt.figure(figsize=(8,12))

    ax1 = plt.subplot2grid((4, 2), (0, 0), colspan=1, rowspan=2)  # 1st row, 1st column, spanning 1 column
    ax2 = plt.subplot2grid((4, 2), (0, 1), colspan=1, rowspan=2)  # 1st row, 2nd column, spanning 1 column
    ax3 = plt.subplot2grid((4, 2), (2, 0), colspan=2)  # 2nd row, 1st column, spanning 2 columns
    ax4 = plt.subplot2grid((4, 2), (3, 0), colspan=2)  # 3rd row, 1st column, spanning 2 columns

    if pattern==0:
        ax1.scatter(electrode_1_zy[0], electrode_1_zy[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)
    elif pattern==4:
        ax1.scatter(electrode_1_zy[0], electrode_1_zy[1], color='gold', s=110, marker='s', label='Return electrode in L4', zorder=3)
    elif pattern==5:
        ax1.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
    elif pattern==6:
        ax1.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
    elif pattern==8:
        ax1.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
        ax1.scatter(electrode_3_zy[0], electrode_3_zy[1], color='yellow', s=110, marker='s', label='Return electrode 2 in L2/3', zorder=3)
    else:
        ax1.scatter(electrode_2_zy[0], electrode_2_zy[1], color='gold', s=110, marker='s', label='Return electrode 1 in L2/3', zorder=3)
        ax1.scatter(electrode_3_zy[0], electrode_3_zy[1], color='yellow', s=110, marker='s', label='Return electrode 2 in L2/3', zorder=3)

        
    ax1.axline(electrode_0_zy, electrode_1_zy, color='limegreen', label='Along layer')
    ax1.axline(electrode_0_zy, electrode_2_zy, color='darkgreen', label='Along column')
    ax1.scatter(electrode_0_zy[0], electrode_0_zy[1], color='orange', s=110, marker='s', label='Central electrode', zorder=3)
    ax1.scatter(node_pos[:,1], node_pos[:,0], s=90, c="blue", alpha=n_spikes_norm)
    ax1.scatter(max_z_axis, electrode_0_zy[1], color='red', marker='*', s=120, label='Max density', zorder=3)
    ax1.scatter(electrode_0_zy[0], max_y_axis, color='red', marker='*', s=120, zorder=3)
    ax1.scatter(max_z_axis,max_y_axis, color='red', marker='*', s=120, zorder=3)

    ax1.set_xlabel('Z Coordinate')
    ax1.set_ylabel('Y Coordinate')
    #ax1.set_xlim([-400,400])
    #ax1.set_ylim([100, 800])
    ax1.set_xlim([-400,400])
    ax1.set_ylim([0, 800])
    ax1.invert_yaxis()  # Invert y-axis for better comparison 
    ax1.invert_xaxis()
    ax1.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    ax1.legend(fontsize='8', loc='center left', bbox_to_anchor=(1, 0.5))

    pcm = ax2.pcolormesh(grid_z_2D, grid_y_2D, density_2D, shading='auto')
    ax2.scatter(node_pos[:,1], node_pos[:,0], c=n_spikes, cmap='viridis', edgecolors='k', linewidths=1)
    fig.colorbar(pcm, ax=ax2, label='Values')
    ax2.set_xlabel('Z Coordinate')
    ax2.set_ylabel('Y Coordinate')
    ax2.set_xlim([-400, 400])
    ax2.set_ylim([0, 800])
    ax2.invert_yaxis()  # Invert y-axis for better comparison
    ax2.invert_xaxis()
    ax2.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    ax2.set_title('2D Kernel Density')

    ax3.plot(grid_z, density_z, color='red', linestyle='-')
    ax3.set_xlabel('Distance along the z axis')
    ax3.set_ylabel('Density')
    ax3.set_title('1D Kernel Density Estimate along z axis')

    ax4.plot(grid_y, density_y, color='red', linestyle='-')
    ax4.set_xlabel('Distance along the y axis')
    ax4.set_ylabel('Density')
    ax4.set_title('1D Kernel Density Estimate along y axis')

    pattern_title="Parallel to cortical columns. Pattern"+str(pattern)+". M"+str(mouse)+". Amplitude "+ str(amplitude)+"."
    fig.suptitle(pattern_title)
    plt.tight_layout(h_pad=4)
    plt.savefig('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_2/plots_column/exp_2_full_kde_p'+str(pattern)+'_amp'+str(amplitude)+'_m_'+str(mouse)+'.png')
    #plt.show()
    return max_y_axis, max_z_axis

def plot1_kde(node_pos, n_spikes, pattern, mouse, amplitude):
    grid_y, grid_z, density_y, density_z = projected_kernel_density_estimate(node_pos, n_spikes)
    max_y_axis=grid_y[np.argmax(density_y)][0]
    max_z_axis=grid_z[np.argmax(density_z)][0]

    node_pos= node_pos[:,1:]
    max_spikes=np.max(n_spikes)
    #print("max number spikes", max_spikes)
    n_spikes_norm=n_spikes/max_spikes
    #print(n_spikes_norm)

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
    plt.scatter(max_z_axis, electrode_0_zy[1], color='red', marker='*', s=120, label='Max density', zorder=3)
    plt.scatter(electrode_0_zy[0], max_y_axis, color='red', marker='*', s=120, zorder=3)
    plt.scatter(max_z_axis,max_y_axis, color='red', marker='*', s=120, zorder=3)

    plt.xlabel('Z Coordinate')
    plt.ylabel('Y Coordinate')
    #plt.set_xlim([-400,400])
    #plt.set_ylim([100, 800])
    plt.xlim([-400, 400])
    plt.ylim([0, 800])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
    plt.legend(fontsize='12', loc='upper right')

    pattern_title="Parallel to cortical columns. Pattern"+str(pattern)+". M"+str(mouse)+". Amplitude "+ str(amplitude)+"."
    plt.title(pattern_title)
    plt.savefig('/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_2/plots_column/exp_2_1d_kde_p'+str(pattern)+'_m_'+str(mouse)+'a_'+str(amplitude)+'.png')
    #plt.close()
    #plt.show()
    return max_y_axis, max_z_axis  

#path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
exp=2
pattern_A=7
mouse_A=2
amplitude_A=10
node_pos_A, n_spikes_A = get_spikes(exp=exp,pattern=pattern_A,mouse=mouse_A,amplitude=amplitude_A)

pattern_B=8
mouse_B=2
amplitude_B=10
node_pos_B, n_spikes_B = get_spikes(exp=exp,pattern=pattern_B,mouse=mouse_B,amplitude=amplitude_B)

positions_filtered_A, spikes_filtered_A, threshold_A = filter_spikes(node_pos_A, n_spikes_A)
positions_filtered_B, spikes_filtered_B, threshold_B = filter_spikes(node_pos_B, n_spikes_B)

for pattern in [0,4,5,6,7,8]:
    pattern_1=pattern
    amplitude_1 = 10
    for mouse in [0,1,2]:
        mouse_1=mouse
        node_pos_1, n_spikes_1 = get_spikes(exp=exp,pattern=pattern_1,mouse=mouse_1,amplitude=amplitude_1)
        positions_filtered_1, spikes_filtered_1, threshold_1 = filter_spikes(node_pos_1, n_spikes_1)
        max_y_axis_1, max_z_axis_1 = plot1_kde(positions_filtered_1, spikes_filtered_1, pattern_1, mouse_1,amplitude_1)
        max_y_1,max_z_1 = full_kde(positions_filtered_1, spikes_filtered_1, pattern_1,mouse_1,amplitude_1)

#spearman, pvalue = Spearmancorr(n_spikes_A= n_spikes_A, n_spikes_B=n_spikes_B, pattern_A=pattern_A, pattern_B=pattern_B, threshold_A = threshold_A, threshold_B = threshold_B)
#statistic, pvalue = Pearsoncorrel(n_spikes_A= n_spikes_A, n_spikes_B=n_spikes_B, pattern_A=pattern_A, pattern_B=pattern_B, threshold_A = threshold_A, threshold_B = threshold_B)
#coordin_A, n_spikes_A, y_grid_A, z_grid_A, density_A = kernel_density_estimate(node_pos=node_pos_A,n_spikes=n_spikes_A, pattern=pattern_A)
#grid_y_A, grid_z_A, density_y_A, density_z_A = projected_kernel_density_estimate(node_pos_A, n_spikes_A)
#max_y_A,max_z_A = full_kde(positions_filtered_A, spikes_filtered_A, pattern_A,mouse_A,amplitude_A)

#Underneath: test_code

#coordinates= node_pos_A[:,1:]

'''non_zero_indices = np.nonzero(n_spikes_A)
non_zero_coordinates = coordinates[non_zero_indices]
non_zero_n_spikes_A = n_spikes_A[non_zero_indices]
# Plot only the non-zero values with color scale based on the values in non_zero_n_spikes_A
plt.scatter(
    non_zero_coordinates[:, 1],  # x-coordinates of the non-zero data points
    non_zero_coordinates[:, 0],  # y-coordinates of the non-zero data points
    c=non_zero_n_spikes_A,       # color scale based on the values in non_zero_n_spikes_A
    cmap='viridis',              # colormap used for coloring the points
    edgecolors='black',          # color of the edges of the points
    linewidths=1                 # width of the edges of the points
)
plt.colorbar(label='n_spikes_A')  # Add colorbar indicating the values of n_spikes_A
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Non-Zero Values from n_spikes_A')
plt.show()'''

'''plt.scatter(
    coordinates[:, 1],      # x-coordinates of the data points
    coordinates[:, 0],      # y-coordinates of the data points
    c=np.where(n_spikes_A == 0, 'white', 'red'),  # conditional color based on n_spikes_A values
    edgecolors='black',     # color of the edges of the points
    linewidths=1            # width of the edges of the points
)
plt.title("plot spiking neurons")
plt.show()'''

#print(p_value_wilcoxon)
#p_value_test = discriminate_signed_rank(n_spikes_A=[0,2,3,0], n_spikes_B=[0,1,3,0], pattern_A=0, pattern_B=4)
#print(p_value_test)
#print(n_spikes_A[149472])
#print(n_spikes_B[149472])
#print(n_spikes_A.shape)
#print(node_pos_A.shape)
#print(n_spikes_B.shape)
#print(n_spikes_A)
#print(n_spikes_B)
#print(node_pos_B.shape)
#print(np.max(n_spikes_A))
#print(np.max(n_spikes_B))
