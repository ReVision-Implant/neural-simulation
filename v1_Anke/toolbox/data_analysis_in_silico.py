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
from scipy.optimize import curve_fit


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

def electrode_coordin(exp,pattern):
    el_coordinates = {
                    0: [-9, 300, 16], # z, y, x coordin
                    1: [-9, 300, 198],
                    2: [-9, 170, 16],
                    3: [-9, 170, 380],
                    4: [-9, 170, 380],
                    5: [-9, 300, 380],
                    6: [-9, 170, -166],
                    7: [-9, 170, -348],
                    8: [-9, 300, -348]
                }

    if exp==4:
        exp_= {
        0: [0,1],
        5: [0,2],
        #7: [0,2,3]
            }
    else:
        exp_= {
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
    
    electrodes= exp_[pattern]
    return_el= electrodes[1]
    location_central_el = np.array([-9, 300, 16])
    #np.transpose(location_central_el)
    location_return= np.array(el_coordinates[return_el])
    #np.transpose(location_return)
    return location_central_el, location_return

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

def project_neurons(exp, pattern, mouse, amplitude, plot=False):
    node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    #print("shape node pos", node_pos.shape)
    pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)
    #print("shape of pos_filtered is ", pos_filtered.shape)
    #print("spikes filtered", spikes_filtered)
    point1, point2 = electrode_coordin(exp, pattern)
    #print("shape point 1", point1.shape)
    #print("point 1", point1)

    point1_2d = np.array([point1[2],point1[1]])
    point2_2d = np.array([point2[2],point2[1]])
    #print("point1_2d", point1_2d)
    #print("point2 2d", point2_2d)
    #print("shape point 1 2d is", point1_2d.shape)

    direction_vector =  np.array([point2_2d[0] - point1_2d[0], point2_2d[1] - point1_2d[1]])

    #print("shape of direction vector is ", direction_vector.shape)
    #print("direction vector", direction_vector)

    pos_filtered_2d = np.column_stack((pos_filtered[:, 2], pos_filtered[:, 1]))
    #print("shape of pos_filtered_2d is ", pos_filtered_2d.shape)
    #print( "pos filtered 2d are", pos_filtered_2d)


    projected_points = np.array([np.dot(np.array(i) - np.array(point1_2d), direction_vector) / np.dot(direction_vector, direction_vector) for i in pos_filtered_2d])

    projected_points, n_spikes_sorted = zip (*sorted(zip(projected_points, spikes_filtered))) # Sort the projected_points in ascending order
    projected_points = np.array(projected_points)
    n_spikes_sorted = np.array(n_spikes_sorted)
    #print("projected points shape is ", projected_points.shape)

    if plot==True:
            fig, ax = plt.subplots(2, 1)
            ax = ax.flatten()

            # Plot for original points, line, and projected points
            ax[0].plot([point1_2d[0], point2_2d[0]], [point1_2d[1], point2_2d[1]], color='blue', label='Line')
            #print( max(n_spikes_sorted))
            for cell in range(pos_filtered_2d.shape[0]):
                projected_coordinate = [point1_2d[0] + projected_points[cell] * direction_vector[0], point1_2d[1] + projected_points[cell] * direction_vector[1]]
                ax[0].scatter(pos_filtered_2d[cell][0], pos_filtered_2d[cell][1], s=50, c="blue", alpha=n_spikes_sorted[cell]/max(n_spikes_sorted))
                ax[0].scatter(projected_coordinate[0], projected_coordinate[1], s=20, c="red", alpha=n_spikes_sorted[cell]/max(n_spikes_sorted))
                #print(pos_filtered_2d[cell][0], pos_filtered_2d[cell][1], n_spikes_sorted[cell]/max(n_spikes_sorted))
            ax[0].set_xlabel('X Coordinate')
            ax[0].set_ylabel('Y Coordinate')
            #ax[0].set_xlim([0, 397])
            #ax[0].set_ylim([0, 380])
            ax[0].invert_yaxis()  # Invert y-axis for better comparison with ImageJ
            ax[0].set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
            ax[0].set_title('Projection of points onto projection axis')
            ax[0].legend()

            # Plot for points along the projection axis
            ax[1].scatter(projected_points, n_spikes_sorted, color='blue')
            ax[1].set_xticks([])
            ax[1].set_ylabel('number of spikes')
            ax[1].set_title('Points along projection axis')
            ax[1].set_aspect('auto')

            fig.tight_layout()
            plt.show()

    return projected_points, n_spikes_sorted

def fit_neurons_kde(projected_points, sorted_spikes, plot=False):
    kde = KernelDensity(bandwidth=200, kernel='gaussian')
    kde.fit(projected_points.reshape(-1,1), sample_weight=sorted_spikes) 

    # Define grid along the projected axis
    grid_size = 100
    grid = np.linspace(min(projected_points), max(projected_points), grid_size).reshape(-1, 1)
    density = np.exp(kde.score_samples(grid))
    centroid = grid[np.argmax(density)][0]

    if plot==True:
            plt.figure()
            # Plot for points along the projection axis
            plt.scatter(projected_points, sorted_spikes, color='blue')
            plt.plot(grid, density, color='red', linestyle='-')
            plt.scatter(centroid, 0, color='red', linewidth=10)
            plt.xlabel('Density along projected axis')
            plt.xticks([0,1],["Point 1", "Point 2"])
            plt.title('1D KDE of projected points')
            plt.show()

    return centroid

def fit_neurons_stdev(projected_points, sorted_spikes, plot=False):
    centroid = np.average(projected_points, weights=sorted_spikes)
    stdev = np.std(projected_points)

    if plot==True:
            plt.figure()
            # Plot for points along the projection axis
            plt.scatter(projected_points, sorted_spikes, color='blue')
            plt.scatter(centroid, 0, color='red')
            plt.plot([centroid-stdev, centroid+stdev], [0, 0], color='red')
            plt.xlabel('Density along projected axis')
            plt.xticks([0,1],["Point 1", "Point 2"])
            plt.title('Centroid and standard deviation of projected points')
            plt.show()
    
    return centroid,stdev

def gauss(x, x0, sigma, amp, offset):
    return offset + amp * np.exp(-(x-x0)**2/(2*sigma**2))

def fit_neurons_gaussian(projected_points, sorted_spikes, plot = False):

    # Take the mean and standard deviation as initial guess for the peak and standard deviation of the Gaussian
    x0_estimated, mu0_estimated = fit_neurons_stdev(projected_points, sorted_spikes, plot=False)

    nb_bins = 20 # or set to a fixed size
    # Create 10 equidistant bins
    bins = np.linspace(projected_points[0], projected_points[-1], nb_bins+1)  # 11 edges for 10 bins
    cumulative_values = np.zeros(nb_bins)

    bin_index = 0
    # Accumulate the spike rates in the bins
    for i in range(len(projected_points)):
        # Check if the projected point exceeds the current bin edge
        while bin_index < nb_bins-1 and projected_points[i] > bins[bin_index + 1]:
            bin_index += 1
        # Add the spike rate value to the corresponding bin
        cumulative_values[bin_index] += sorted_spikes[i]

    # Prepare the new list with mid-points of the bins and their corresponding cumulative values
    points = (bins[:-1] + bins[1:]) / 2
    amplitude = cumulative_values

    p0 = [x0_estimated, mu0_estimated, max(amplitude), min(amplitude)] # Initial guess for p = [mu, sigma, amplitude, offset]
    p_bounds = ((x0_estimated*0.8 if x0_estimated >= 0 else x0_estimated*1.2, 0, min(amplitude), min(amplitude)*0.5), (x0_estimated*1.2 if x0_estimated >= 0 else x0_estimated*0.8, np.inf, max(amplitude)*2, max(amplitude))) # Boundaries for the different parameters
    popt, pcov = curve_fit(gauss, points, amplitude, p0 =p0, bounds=p_bounds, maxfev=5000, method='trf', loss='linear') # Fit the points and their values to the given function todo possibly choose different loss function for regularization against outliers (eg cauchy)
    fitted_curve = gauss(points, popt[0], popt[1], popt[2], popt[3]) # Compute the fitted curve for the projected_points

    if plot==True:
        plt.figure()
        # Plot for points along the projection axis
        plt.scatter(points, amplitude, color='blue')
        plt.plot(points, fitted_curve, color='red', linestyle='-')
        plt.xlabel('Density along projected axis')
        plt.xticks([0,1],["Point 1", "Point 2"])
        plt.title('1D Gaussian fit of projected points')
        plt.show()

    return popt

def spatial_analysis(exp, pattern, mouse, amplitude):
    node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)

    # Centroid along cortical column:
    projected_coordinates, n_spikes_sorted = project_neurons(exp, pattern, mouse, amplitude)
    centroid_along_column = fit_neurons_kde(projected_coordinates, n_spikes_sorted, plot=False)

    # Centroid along cortical layer
    [centroid_along_layer, stdev_along_layer, gaussian_amp, gaussian_offset] = fit_neurons_gaussian(projected_coordinates, n_spikes_sorted, plot=plot)
    
    return centroid_along_column, centroid_along_layer, stdev_along_layer


################################################## TEST CODE ##########################################################
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

#correlation_per_angle(exp=[4,4], patterns=[0,5], mouse=[0,0], amplitude=[10,10])

projected_neurons,spikes_sorted = project_neurons(exp=4, pattern=0, mouse=0, amplitude=10, plot = False)
#centroid = fit_neurons_kde(projected_neurons, spikes_sorted, plot=False)
#centroid, stdev = fit_neurons_stdev(projected_neurons, spikes_sorted, plot = True)
#popt_test = fit_neurons_gaussian(projected_neurons, spikes_sorted, plot = True)


############################################################
######### DATA ANALYSIS DIRECTIONALITY SPATIAL #############
############################################################