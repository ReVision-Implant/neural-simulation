import pandas as pd
from hdf5 import HDF5
import numpy as np
import os
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

def get_spikes(exp,pattern,mouse,amplitude, v1=True, **kwargs):
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

    path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
    

    nodes_dirs= [str(path)+'/virtual_mice_mask/mouse_'+str(mouse)+'/v1_nodes.h5']
    spikes_dirs= [str(path)+'/exp_'+str(exp)+'/output/pattern_'+str(pattern)+'/amplitude_'+str(amplitude)+'/mouse_'+str(mouse)+'/spikes.csv']
    spikes_bkg_dirs= [str(path)+'/exp_'+str(exp)+'/output/bkg/mouse_'+str(mouse)+'/spikes.csv']
        
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

def discriminate_signed_rank(n_spikes_A, n_spikes_B,pattern_A,pattern_B):
        '''
        Use the Wilcoxon signed-rank test to get a p-value as index of separability between the two neuronal populations.
        '''
        signed_rank = wilcoxon(n_spikes_A, n_spikes_B) # Apply Wilcoxon test
        print('P-value for Wilcoxon signed-rank test for stim patterns '+str(pattern_A)+' and '+str(pattern_B)+' is ' + str(round(signed_rank.pvalue,5)))
        return(signed_rank.pvalue)

def kernel_density_estimate(node_pos, n_spikes, pattern):
        '''
        2D Kernel Density Estimate of the data
        '''
        non_zero_indices = np.nonzero(n_spikes)
        coordinates= node_pos[:,1:] #select only the y and z coordinates
        coordinates = coordinates[non_zero_indices]
        n_spikes= n_spikes[non_zero_indices]

        kde = KernelDensity(bandwidth=10, kernel='gaussian') # Choose model and parameters
        ###vanaf hier verder werken
        kde.fit(coordinates, sample_weight=n_spikes) # Train model

        grid_size = 100 # 100 points in x and in y direction
        y_grid, z_grid = np.meshgrid(np.linspace(100, 800, grid_size), np.linspace(-250, 500, grid_size))
        grid_points = np.vstack([y_grid.ravel(), z_grid.ravel()]).T
        density = np.exp(kde.score_samples(grid_points)).reshape(y_grid.shape) # Evaluate model for all points on the grid

        fig = plt.figure()
        plt.pcolormesh(z_grid, y_grid, density, shading='auto')
        plt.scatter(coordinates[:,1], coordinates[:,0], c=n_spikes, cmap='viridis', edgecolors='k', linewidths=1)
        plt.colorbar(label='Values')
        plt.xlabel('Z Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([-250, 400])
        plt.ylim([100, 800])
        #plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        # plt.legend()
        plt.title('Kernel Density Estimate for stim. pattern ' + str(pattern))
        plt.show()
        #plt.close()

        return coordinates, n_spikes, y_grid, z_grid, density

def projected_kernel_density_estimate(node_pos,n_spikes):
        '''
        Projection of the data before the Kernel Density Estimate is useful when trying the understand the spatial
        distribution of cell activity along a specific axis. The density estimate now reflects the distribution
        of the data along the projected data.
        If projection would only take place after the Kernel Density Estimate, then you just integrate the density
        function, but the density estimate is still based on a 2D-distribution.
        '''
        non_zero_indices = np.nonzero(n_spikes)
        coordinates= node_pos[:,1:] #select only the y and z coordinates
        coordinates = coordinates[non_zero_indices]
        n_spikes= n_spikes[non_zero_indices]

        projected_points_z=coordinates[:,1]
        projected_points_y=coordinates[:,0]
        # Define grid along the projected axis
        grid_size = 100
        grid_z = np.linspace(min(projected_points_z), max(projected_points_z), grid_size).reshape(-1, 1)
        grid_y = np.linspace(min(projected_points_y), max(projected_points_y), grid_size).reshape(-1, 1)

        # Perform kernel density estimation
        kde = KernelDensity(bandwidth=5, kernel='gaussian') 

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
        plt.show()
        #plt.close()

        return grid_y, grid_z, density_y, density_z

path ='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke'
exp=2
pattern_A=0
mouse=1
amplitude=10
node_pos_A, n_spikes_A = get_spikes(exp=exp,pattern=pattern_A,mouse=mouse,amplitude=amplitude)

pattern_B=4
node_pos_B, n_spikes_B = get_spikes(exp=exp,pattern=pattern_B,mouse=mouse,amplitude=amplitude)

#p_value_wilcoxon = discriminate_signed_rank(n_spikes_A= n_spikes_A, n_spikes_B=n_spikes_B, pattern_A=0, pattern_B=4)
#coordin_A, n_spikes_A, y_grid_A, z_grid_A, density_A = kernel_density_estimate(node_pos=node_pos_A,n_spikes=n_spikes_A, pattern=pattern_A)
grid_y_A, grid_z_A, density_y_A, density_z_A = projected_kernel_density_estimate(node_pos_A, n_spikes_A)

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

#p_value_test = discriminate_signed_rank(n_spikes_A=[0,2,3,0], n_spikes_B=[0,1,3,0])
#print(p_value_test)
#print(p_value_wilcoxon)
#print(n_spikes_A[149472])
#print(n_spikes_B[149472])
#print(n_spikes_A.shape)
#print(node_pos_A.shape)
#print(n_spikes_B.shape)
#print(node_pos_B.shape)
#print(np.max(n_spikes_A))
#print(np.max(n_spikes_B))
