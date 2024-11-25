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
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.ticker import FuncFormatter
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

        spikes = pd.read_csv(spikes_dir, sep=r'\s+')
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
        
        y_coordin=node_pos[:,1]
        n_spikes_L234=[]
        node_pos_L234=[]
        for index, y in enumerate(y_coordin):
            if y >= 100 and y<=430: # only select neurons in layer 2/3 and 4 of the cortex
                n_spikes_L234.append(n_spikes[index])
                node_pos_L234.append(node_pos[index,:])
        n_spikes_L234=np.array(n_spikes_L234)
        node_pos_L234=np.array(node_pos_L234)
        #print(n_spikes_L234.shape)
        #print(node_pos_L234.shape)
    
    #for analysis all neurons:
    return node_pos, n_spikes

    #only looking at neurons in layer 2/3 and 4
    #return node_pos_L234, n_spikes_L234

#test code
#n_spikes_L234_test, node_pos_L234_test = get_spikes(4,0,0,10)

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
        if value >= threshold:
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
    
        #print('The Spearman correlation coefficient for stim patterns ' + str(pattern_A) + ' and ' + str(pattern_B) +
        #  ' is ' + str(round(statistic, 2)) + ', the p-value is ' + str(round(pvalue, 4)) +
        #  ', and the 95% confidence interval is [' + str(round(lower_bound, 2)) + ', ' + str(round(upper_bound, 2)) + ']')

        return(statistic)

def overlap(n_spikes_A, n_spikes_B, threshold_A, threshold_B): 
        
        #print("threshold pattern 1", threshold_A, " threshold pattern 2:", threshold_B)
        #n_spikes_A_filtered=[]
        #n_spikes_B_filtered=[]
        activity_A =[]
        activity_B =[]
        for value1, value2 in zip(n_spikes_A, n_spikes_B):
            if value1 >= threshold_A or value2 >=threshold_B:
            #print(value1,value2)
                #n_spikes_A_filtered.append(value1)
                #n_spikes_B_filtered.append(value2)

                if value1 >= threshold_A:
                    activity_A.append(1)
                else:
                    activity_A.append(0)
                if value2 >= threshold_B:
                    activity_B.append(1)
                else:
                    activity_B.append(0)     

        #n_spikes_A= n_spikes_A_filtered
        #n_spikes_B= n_spikes_B_filtered

        #print(threshold_A)
        #print(n_spikes_A)
        #print(activity_A)

        overlap = [1 if (activity_A[i] == 1 and activity_B[i] == 1) else 0 for i in range(len(activity_A))]

        active_neurons_total = len(activity_A)
        #active_neurons_A = activity_A.count(1)
        #active_neurons_B = activity_B.count(1)
        active_neurons_overlap = overlap.count(1)
        #print("active neurons overlap", active_neurons_overlap)
        #print("total active neurons", active_neurons_total)

        if active_neurons_total > 0:
            return np.round(active_neurons_overlap/active_neurons_total,4), active_neurons_total
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
                    3: [-9, 170, 198],
                    4: [-9, 170, 380],
                    5: [-9, 300, 380],
                    6: [-9, 170, -166],
                    7: [-9, 170, -348],
                    8: [-9, 300, -166],
                    9: [-9, 300, -348],
                    10: [-9, 170, 289 ], #in between el 3 and 4
                    11: [-9, 170, 107], #in between el 2 and 3
                    12: [-9, 170, -75], #in between el 2 and 6
                    13: [-9, 170, -257] #in betzeen el 6 and 7

                }

    if exp==4:
        exp_= {
        0: [0,1],
        5: [0,2],
        #7: [0,2,3]
        7: [0,11]
            }
    else:
        exp_= {
            0: [0,5],
            1: [0,4],
            #2: [0,3,4],
            2: [0,10],
            3: [0,6],
            4: [0,7],
            5: [0,8],
            6: [0,9],
            #7: [0,2,6],
            7: [0,12],
            #8: [0,6,7],
            8: [0,13],
            9: [0,3]
                    } 
    
    electrodes= exp_[pattern]
    return_el= electrodes[1]
    location_central_el = np.array([-9, 300, 16])
    #np.transpose(location_central_el)
    location_return= np.array(el_coordinates[return_el])
    #np.transpose(location_return)
    return location_central_el, location_return

def correlation_per_angle(exp=[4,5], patterns=[0,1], mouse=[0,0], amplitude=10):
        angles = []
        correlations = []
        overlaps = []
    
        for i in range(len(patterns)):
            node_pos_i, n_spikes_i= get_spikes(exp[i],patterns[i],mouse[i],amplitude)
            pos_filtered_i, spikes_filtered_i, threshold_i = filter_spikes(node_pos_i,n_spikes_i)
            
            for j in range(i+1,len(patterns)):
                node_pos_j, n_spikes_j= get_spikes(exp[j],patterns[j],mouse[j],amplitude)
                pos_filtered_j, spikes_filtered_j, threshold_j = filter_spikes(node_pos_j,n_spikes_j)

                #correlations
                statistic = correlation(n_spikes_i, n_spikes_j,patterns[i],patterns[j],threshold_i, threshold_j)
                correlations.append(statistic)

                #overlaps
                overlap_= overlap(n_spikes_i, n_spikes_j, threshold_i, threshold_j)
                overlaps.append(overlap_)

                #angles
                location_central_el, location_return_i = electrode_coordin(exp[i], patterns[i])
                location_central_el, location_return_j = electrode_coordin(exp[j], patterns[j])
                angle= get_electrode_angles(location_central_el, location_return_i,location_return_j)
                angles.append(angle)
                #print('angle between exp', exp[i], 'p', patterns[i], 'and exp', exp[j], 'p', patterns[j], 'is', angle, "with an overlap ", overlap_[0])
        angles, correlations, overlaps = zip(*sorted(zip(angles, correlations, overlaps))) # Sort the correlations in ascending order
        return angles, correlations, overlaps

def project_neurons(exp, pattern, mouse, amplitude, point1= [-9, 300, 16], point2=[-9, 170, 16], plot=False):
    node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    #print("shape node pos", node_pos.shape)
    pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)
    #print("shape of pos_filtered is ", pos_filtered.shape)
    #print("spikes filtered", spikes_filtered)

    #when projecting on the stimulation direction
    #point1, point2 = electrode_coordin(exp, pattern)

    #print("shape point 1", point1.shape)
    #print("point 2", point2)

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

    #normalisation
    projected_points = np.array([np.dot(np.array(i) - np.array(point1_2d), direction_vector) / np.dot(direction_vector, direction_vector) for i in pos_filtered_2d])
    #print("projected_points", projected_points)
    #no normalisation
    #projected_points= np.array([np.dot(np.array(i) - np.array(point1_2d), direction_vector) for i in pos_filtered_2d])


    projected_points, n_spikes_sorted = zip (*sorted(zip(projected_points, spikes_filtered))) # Sort the projected_points in ascending order
    projected_points = np.array(projected_points)
    n_spikes_sorted = np.array(n_spikes_sorted)
    #print("projected_points", projected_points)
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

    #print(popt) 

    return popt

def spatial_analysis(exp, pattern, mouse, amplitude):
    #node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    #pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)

    # Centroid along cortical column:
    point1 = [-9, 300, 16]
    point2 = [-9, 170, 16]
    projected_coordinates, n_spikes_sorted = project_neurons(exp, pattern, mouse, amplitude, point1, point2)
    centroid_along_column = fit_neurons_kde(projected_coordinates, n_spikes_sorted, plot=False)

    # Centroid along cortical layer
    point3= [-9, 300, 198]
    projected_coordinates, n_spikes_sorted = project_neurons(exp, pattern, mouse, amplitude, point1, point3)
    [centroid_along_layer, stdev_along_layer, gaussian_amp, gaussian_offset] = fit_neurons_gaussian(projected_coordinates, n_spikes_sorted, plot=False)
    
    #print("centroid_along_column", centroid_along_column, "centroid_along_layer", centroid_along_layer)
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

#projected_neurons,spikes_sorted = project_neurons(exp=4, pattern=0, mouse=0, amplitude=20, plot = False)
#centroid = fit_neurons_kde(projected_neurons, spikes_sorted, plot=False)
#centroid, stdev = fit_neurons_stdev(projected_neurons, spikes_sorted, plot = True)
#popt_test = fit_neurons_gaussian(projected_neurons, spikes_sorted, plot = False)

#centroid_c, centroid_l, stdev = spatial_analysis(exp=4,pattern=0,mouse=0, amplitude=10)


############################################################
######### DATA ANALYSIS DIRECTIONALITY SPATIAL #############
############################################################

def electrode_dir(exp,pattern):
    if exp==4:
        electrode_dir= {
        0: (-1,0),
        5: (0,1),
        7: (-0.5,1)
            }
    else:
        electrode_dir= {
            0: (-2,0),
            1: (-2,1),
            2: (-1.5,1),
            3: (1,1),
            4: (2,1),
            5: (1,0),
            6: (2,0),
            7: (0.5,1),
            8: (1.5,1),
            9: (-1,1)
                    } 
    
    electrodes_direction= electrode_dir[pattern]
    return electrodes_direction

def directionality_spatial(exp=[4,5], patterns=[0,0], mice=[0,1,2], amplitude=10):
    #print("len exp",len(exp))
    colors = ['red', 'blue', 'orange', 'purple', 'yellow', 'grey', 'green', 'pink', 'brown', 'cyan', 'magenta', 'teal', 'lime']
    centroids = [[] for _ in range(len(exp) * len(mice))]
    for pattern in range (0,len(exp)):
          for i in range(0,len(mice)):
            #print(exp[pattern])
            #print("exp", exp[pattern], "pattern", patterns[pattern], "mouse", mice[i])
            centroid_column, centroid_layer, stdev = spatial_analysis(exp[pattern], patterns[pattern], mice[i], amplitude)
           # print("centroid column", centroid_column)
            #print("centroid layer", centroid_layer)
            centroids[pattern * len(mice) + i].append([centroid_column, centroid_layer, stdev])
    #print(centroids)  

    plt.figure
    for pattern in range(0,len(exp)):
        for i in range(0, len(mice)):
            #centroid_collection = []
            for centroid in range(len(centroids[i])):
                plt.plot(centroids[pattern*len(mice)+i][centroid][1], centroids[pattern*len(mice)+i][centroid][0], 'o', markersize=5, color=colors[pattern], alpha=0.5, label='Intra-slice centroids')
            
                #print(str(centroids[pattern*len(mice)+i][centroid][1]), colors[pattern])
                # Plot a polygon for illustration purposes
                #centroid_collection.append([centroids[i][centroid][1], centroids[i][centroid][0]])

        plt.plot(np.mean([centroids[pattern*len(mice) +j][centroid][1] for j in range(len(mice))]), 
                 np.mean([centroids[pattern*len(mice) +j][centroid][0] for j in range(len(mice))]), 'o', markersize=10, color=colors[pattern], label='Inter-slice centroid', zorder=3)
        
        return_elec = electrode_dir(exp[pattern], patterns[pattern])
        plt.arrow(0, 0, return_elec[0]*-1, return_elec[1], length_includes_head=True, width=0.01, head_width=0.03, color=colors[pattern])
    
        # Plot a polygon for illustration purposes
        #hull = ConvexHull(np.array(centroid_collection))
        #for simplex in hull.simplices:
        #    plt.plot([centroid_collection[simplex[0]][0], centroid_collection[simplex[1]][0]], [centroid_collection[simplex[0]][1], centroid_collection[simplex[1]][1]], color=colors[i], alpha=0.6)
        #plt.fill([centroid_collection[vertix][0] for vertix in hull.vertices], [centroid_collection[vertix][1] for vertix in hull.vertices], color=colors[i], alpha=0.2)

    plt.xlabel('Distance parallel to surface [µm]', fontsize=20)
    plt.ylabel('Distance perpendicular \nto surface [µm]', fontsize=20)
    plt.gca().invert_xaxis()
    #plt.gca().invert_yaxis()

    plt.xticks(ticks=[-1, -0.5, 0, 0.5, 1], labels=['182', '91', '0', '-91', '-182'], fontsize=15)
    plt.yticks(ticks=[1, 0.66, 0.33, 0, -0.33, -0.66, -1, -1.33, -1.66], labels=['150','200','250','300', '350', '400', '450', '500', '550'], fontsize=15)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
    plt.legend(by_label.values(), by_label.keys(), fontsize=12)
    plt.gca().get_legend().legend_handles[0].set_color('black')
    plt.gca().get_legend().legend_handles[1].set_color('black')

    plt.text(-0.8, 0.2, 'Layer 2/3', fontsize=15)
    plt.hlines(0.07, -0.95, -0.25, colors='black', linestyles='dashed')
    plt.text(-0.8, -0.15, 'Layer 4', fontsize=15)

    plt.title('Activity centroids for different directions', fontsize=20)
    plt.tight_layout()
    plt.show()

#directionality_spatial(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mice=[0,1,2], amplitude=10)

############################################################
######### DATA ANALYSIS ASYMMETRY #############
############################################################

def asymmetry_correlation(mice=[0,1,2]):
    correlatation_asymm_10 = []
    overlap_asymm_10 = []
    
    for i in mice:
        node_pos_A, n_spikes_A = get_spikes(exp=4, pattern=0, mouse=i, amplitude=10)
        positions_filtered_A, spikes_filtered_A, threshold_A = filter_spikes(node_pos_A, n_spikes_A)
        
        node_pos_B, n_spikes_B = get_spikes(exp=4, pattern=4, mouse=i, amplitude=10)
        positions_filtered_B, spikes_filtered_B, threshold_B = filter_spikes(node_pos_B, n_spikes_B)

        correlation_i = correlation(n_spikes_A, n_spikes_B, pattern_A=0, pattern_B=4, threshold_A=threshold_A, threshold_B=threshold_B)
        overlap_i = overlap(n_spikes_A, n_spikes_B, threshold_A, threshold_B)

        correlatation_asymm_10.append(correlation_i)
        overlap_asymm_10.append(overlap_i[0])
        #print(overlap_asymm_10)

    fig, ax = plt.subplots()  

    for correl in correlatation_asymm_10:
        plt.plot(1, correl, 'o', markersize=5, color='orange')
    
    plt.plot(1, np.mean(correlatation_asymm_10), 'o', markersize=10, color='orange')
    ax.errorbar(1, np.mean(correlatation_asymm_10), yerr=np.std(correlatation_asymm_10), capsize=20, capthick=2, color='orange')
    
    plt.ylim([0, 1])
    plt.yticks(fontsize=15)
    plt.ylabel('Spearman correlation coeff.', fontsize=20)
    plt.title('Correlation for symmetry versus asymmetry', fontsize=20)
    fig.tight_layout()

    plt.show()


    # overlap figure does not make too much sense in in-silico exp??
    fig, ax = plt.subplots()  

    for i in overlap_asymm_10:
        plt.plot(1, i, 'o', markersize=5, color='orange')
        #print(i)
    
    plt.plot(1, np.mean(overlap_asymm_10), 'o', markersize=10, color='orange')
    ax.errorbar(1, np.mean(overlap_asymm_10), yerr=np.std(overlap_asymm_10), capsize=20, capthick=2, color='orange')
    #print(np.mean(overlap_asymm_10))
    #print(np.std(overlap_asymm_10))

    plt.ylim([0, 1])
    plt.yticks(fontsize=15)
    plt.ylabel('Binary overlap', fontsize=20)
    plt.title('Overlap for symmetry versus asymmetry', fontsize=20)
    fig.tight_layout()

    #plt.show()

#asymmetry_correlation(mice=[0,1,2])

############################################################
######### DATA ANALYSIS DIRECTIONALITY CELLULAR ############
############################################################

def directionality_cellular_corr(exp=[4,5], patterns=[0,0], mouse=0):
    mouse_x= [mouse]*len(exp)
    markersize = 7
    colors = ['blue', 'orange']
    plt.figure()
    angles_20, correlations_20, overlaps_20 = correlation_per_angle(exp=exp, patterns=patterns, mouse=mouse_x, amplitude=20)
    for combination in range(len(angles_20)):
         plt.plot(angles_20[combination], correlations_20[combination],'o', markersize=markersize, color=colors[0], label= '20 µA')

    angles_10, correlations_10, overlaps_10 = correlation_per_angle(exp=exp, patterns=patterns, mouse= mouse_x, amplitude=10)
    for combination in range(len(angles_10)):
       plt.plot(angles_10[combination], correlations_10[combination],'o', markersize=markersize, color=colors[1], label= '10 µA')  
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
    plt.legend(by_label.values(), by_label.keys(), fontsize=12)

    plt.xlabel('Angle between 2 directions [°]', fontsize=20)
    plt.ylabel('Spearman correlation coeff. \nof 2 directions', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    bottom, top = plt.ylim()
    plt.ylim(bottom, 1)
    plt.title('Correlation for different directions, mouse '+ str(mouse), fontsize=20)
    plt.tight_layout()
    plt.show()

#directionality_cellular_corr(exp=[4,4], patterns=[0,5], mouse=0)
#directionality_cellular_corr(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=0)
#directionality_cellular_corr(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=1)
#directionality_cellular_corr(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=2)

def directionality_cellular_overlap(exp=[4,5], patterns=[0,9], mouse=0):
    import numpy as np
    import matplotlib.pyplot as plt

    mouse_x = [mouse] * len(exp)
    markersize = 7
    colors = ['blue', 'orange']
    plt.figure()
    
    # Generate sample data
    angles_20, correlations_20, overlaps_20 = correlation_per_angle(exp=exp, patterns=patterns, mouse=mouse_x, amplitude=20)
    angles_10, correlations_10, overlaps_10 = correlation_per_angle(exp=exp, patterns=patterns, mouse=mouse_x, amplitude=10)

    # Plotting the data points
    for combination in range(len(angles_20)):
        plt.plot(angles_20[combination], overlaps_20[combination][0], 'o', markersize=markersize, color=colors[0], label='20 µA')
        #print("angles 20", angles_20)
    for combination in range(len(angles_10)):
        plt.plot(angles_10[combination], overlaps_10[combination][0], 'o', markersize=markersize, color=colors[1], label='10 µA')  
        #print("angles 10", angles_10)

    # Filter and sort data for linear fitting
    filtered_sorted_pairs20 = sorted((x, y[0]) for x, y in zip(angles_20, overlaps_20) if not np.isnan(x) and not np.isnan(y[0]))
    filtered_sorted_pairs10 = sorted((x, y[0]) for x, y in zip(angles_10, overlaps_10) if not np.isnan(x) and not np.isnan(y[0]))

    sorted_angles20, sorted_overlaps20 = zip(*filtered_sorted_pairs20)
    sorted_angles10, sorted_overlaps10 = zip(*filtered_sorted_pairs10)

        
    slope20, intercept20 = np.polyfit(sorted_angles20, sorted_overlaps20, 1)
    slope10, intercept10 = np.polyfit(sorted_angles10, sorted_overlaps10, 1)

    # Generate line for plotting
    x_line = np.linspace(0, 180, 100)
    y_line20 = slope20 * x_line + intercept20
    y_line10 = slope10 * x_line + intercept10

    linewidth = 3
    plt.plot(x_line, y_line20, color='blue', linewidth=linewidth)
    equation_text = f"y = {slope20:.5f} x + {intercept20:.2f}"
    plt.text(0.05 * max(x_line), 0.9 * max(y_line20), equation_text, color='black', fontsize=10)

    plt.plot(x_line, y_line10, color='orange', linewidth=linewidth)
    equation_text = f"y = {slope10:.5f} x + {intercept10:.2f}"
    plt.text(0.05 * max(x_line), 0.9 * max(y_line10), equation_text, color='black', fontsize=10)

    # Legend and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Remove duplicate labels
    plt.legend(by_label.values(), by_label.keys(), fontsize=12)

    plt.xlabel('Angle between 2 directions [°]', fontsize=20)
    plt.ylabel('Binary overlap of 2 directions', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    bottom, top = plt.ylim()
    plt.ylim(bottom, 1)
    plt.title('Neural overlap for different directions, mouse ' + str(mouse), fontsize=20)
    plt.tight_layout()
    plt.show()

#directionality_cellular_overlap(exp=[4,4], patterns=[5,7], mouse=0)
#directionality_cellular_overlap(exp=[4,4,4,5,5,5,5], patterns=[0,5,7,9,2,1,0], mouse=0)
#directionality_cellular_overlap(exp=[4,5,5,5,5,5,5], patterns=[5,7,3,8,4,6,5], mouse=0)
#directionality_cellular_overlap(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=0)
#directionality_cellular_overlap(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=1)
#directionality_cellular_overlap(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mouse=2)

def directionality_cellular_overlap_all(exp, patterns, mice):
    import numpy as np
    import matplotlib.pyplot as plt

    markersize = 7
    colors = ['blue', 'orange']
    plt.figure()

    # Data storage for each amplitude
    all_angles_20, all_overlaps_20 = [], []
    all_angles_10, all_overlaps_10 = [], []

    # Loop through each mouse and collect data
    for idx, mouse in enumerate(mice):
        mouse_x = [mouse] * len(exp)
        
        # Generate sample data
        angles_20, correlations_20, overlaps_20 = correlation_per_angle(exp=exp, patterns=patterns, mouse=mouse_x, amplitude=20)
        angles_10, correlations_10, overlaps_10 = correlation_per_angle(exp=exp, patterns=patterns, mouse=mouse_x, amplitude=10)

        # Collect data for amplitude 20
        all_angles_20.extend(angles_20)
        all_overlaps_20.extend(overlap[0] for overlap in overlaps_20)

        # Collect data for amplitude 10
        all_angles_10.extend(angles_10)
        all_overlaps_10.extend(overlap[0] for overlap in overlaps_10)

        # Plotting the data points
        for combination in range(len(angles_20)):
            plt.plot(angles_20[combination], overlaps_20[combination][0], 'o', markersize=markersize, color=colors[0], label='20 µA')
            #print("angles 20", angles_20)
        for combination in range(len(angles_10)):
            plt.plot(angles_10[combination], overlaps_10[combination][0], 'o', markersize=markersize, color=colors[1], label='10 µA')  
            #print("angles 10", angles_10)

    # Filter and sort combined data for amplitude 20
    filtered_sorted_pairs_20 = sorted((x, y) for x, y in zip(all_angles_20, all_overlaps_20) if not np.isnan(x) and not np.isnan(y))
    sorted_angles_20, sorted_overlaps_20 = zip(*filtered_sorted_pairs_20)

    # Filter and sort combined data for amplitude 10
    filtered_sorted_pairs_10 = sorted((x, y) for x, y in zip(all_angles_10, all_overlaps_10) if not np.isnan(x) and not np.isnan(y))
    sorted_angles_10, sorted_overlaps_10 = zip(*filtered_sorted_pairs_10)

    # Perform linear fits
    slope_20, intercept_20 = np.polyfit(sorted_angles_20, sorted_overlaps_20, 1)
    slope_10, intercept_10 = np.polyfit(sorted_angles_10, sorted_overlaps_10, 1)

    # Generate lines for plotting
    x_line = np.linspace(0, 180, 100)
    y_line_20 = slope_20 * x_line + intercept_20
    y_line_10 = slope_10 * x_line + intercept_10

    linewidth = 3
    # Plot linear fits
    plt.plot(x_line, y_line_20, color='blue', linewidth=linewidth)
    plt.plot(x_line, y_line_10, color='orange', linewidth=linewidth)

    # Add equations for fits
    equation_text_20 = f"20 µA: y = {slope_20:.5f}x + {intercept_20:.2f}"
    equation_text_10 = f"10 µA: y = {slope_10:.5f}x + {intercept_10:.2f}"
    plt.text(10, 0.8, equation_text_20, color='black', fontsize=10)
    plt.text(10, 0.7, equation_text_10, color='black', fontsize=10)

    # Legend and labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Remove duplicate labels
    plt.legend(by_label.values(), by_label.keys(), fontsize=12)

    plt.xlabel('Angle between 2 directions [°]', fontsize=20)
    plt.ylabel('Binary overlap of 2 directions', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    bottom, top = plt.ylim()
    plt.ylim(bottom, 1)
    plt.title('Neural overlap for different directions', fontsize=20)
    plt.tight_layout()
    plt.show()

# Call the function for all mice
directionality_cellular_overlap_all(
    exp=[4,4,4,5,5,5,5,5,5,5,5,5,5],
    patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9],
    mice=[0, 1, 2]
)



############################################################
########################## PCA #############################
############################################################

def PCA_analysis(exp=[4, 5], patterns=[0, 0], mice=[0,1,2], amp=10):
    colors = ['blue', 'green', 'red']  

    plt.figure()
    
    for m in mice:
        node_pos_0, n_spikes_0 = get_spikes(exp=exp[0], pattern=patterns[0], mouse=m, amplitude=amp)
        #print("node pos shape", np.array(node_pos_0).shape)
        active_cells = np.zeros(len(node_pos_0))
        #print("active cells shape", np.array(active_cells).shape)
        
        # Mark active cells across all experiments 
        for i in range(len(exp)):
            node_pos_i, n_spikes_i = get_spikes(exp=exp[i], pattern=patterns[i], mouse=m, amplitude=amp)
            positions_filt_i, spikes_filt_i, threshold_i = filter_spikes(node_pos_i, n_spikes_i)
            for cell in range(len(node_pos_i)):
                if n_spikes_i[cell] >= threshold_i:
                    active_cells[cell] += 1

        # Create the activity matrix 
        activity_matrix = []
        for i in range(len(exp)):
            node_pos_i, n_spikes_i = get_spikes(exp=exp[i], pattern=patterns[i], mouse=m, amplitude=amp)
            row = []
            for cell in range(len(n_spikes_i)):
                if active_cells[cell] != 0:
                    row.append(n_spikes_i[cell])
            activity_matrix.append(row)
        
        activity_matrix = np.array(activity_matrix).T  # Transpose to shape (num_active_cells, num_exp)
        #print("shape activity matrix", activity_matrix.shape)

        # Perform PCA 
        pca = PCA(n_components=len(exp))
        pca.fit(activity_matrix)
        explained_variance = pca.explained_variance_ratio_
        
        #print(f"explained var for mouse {m+1}:", explained_variance)
        cumulative_variance = np.cumsum(explained_variance)

        # Plot cumulative variance 
        plt.plot(cumulative_variance, color=colors[m], linewidth=2.5, label=f'Mouse {m + 1}')

    # Add the 90% and 95% variance lines
    plt.hlines(0.95, 0, len(exp) - 1, linestyles='--', colors='gray')
    plt.hlines(0.90, 0, len(exp) - 1, linestyles='--', colors='gray')

    # Add labels, title, and legend
    plt.xlabel('Number of current directions', fontsize=20)
    plt.ylabel('Cumulative variance', fontsize=20)
    plt.gca().get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x+1), ',')))  # Add 1 to every xticklabel
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title('PCA on neural activity for ' + str(len(exp)) + ' directions', fontsize=20)
    plt.legend(fontsize=15)  
    plt.tight_layout()

    plt.show()


#PCA_analysis(exp=[4,4,4,5,5,5,5,5,5,5,5,5,5], patterns=[0,5,7,0,1,2,3,4,5,6,7,8,9], mice=[0,1,2], amp=10)

#########################################################################################
################################### DEPTH ###############################################
#########################################################################################

def point_line_distance(point, line_start, line_end):
        # Line equation: Ax + By + C = 0
        A = line_end[1] - line_start[1]
        B = line_start[0] - line_end[0]
        C = line_end[0] * line_start[1] - line_start[0] * line_end[1]

        # Perpendicular distance formula
        distance = abs(A * point[0] + B * point[1] + C) / math.sqrt(A**2 + B**2)
        return distance

def neurons_on_imaging_plane(exp, pattern, mouse, amplitude, point1= [-9, 300, 16], angle_degrees=60):
    node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    #print("shape node pos", node_pos.shape)
    pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)
    #print("shape of pos_filtered is ", pos_filtered.shape)
    #print("spikes filtered", spikes_filtered)

    angle_radians= math.radians(angle_degrees)
    tangent =math.tan(angle_radians)
    print("tangent",tangent)
    x=(tangent*416)-9

    point2=[416,x]

    #Convert points to 2D (z, x) plane
    point1_2d = np.array([point1[2],point1[0]])
    point2_2d = np.array([point2[1],point2[0]])
    print("point1_2d", point1_2d)
    print("point2 2d", point2_2d)
    #print("shape point 1 2d is", point1_2d.shape)


    y_min = 100
    y_max = 430
    number_neurons_on_implane=0

    distance_threshold = 15.0

    for neuron in range(len(pos_filtered)):
        # Get x, y, z coordinates of the neuron
        x_neuron = pos_filtered[neuron, 0]
        y_neuron = pos_filtered[neuron, 1]
        z_neuron = pos_filtered[neuron, 2]

        # Check if y is within the range
        if y_min <= y_neuron <= y_max:
            neuron_2d = np.array([z_neuron, x_neuron])
            distance_to_line = point_line_distance(neuron_2d, point1_2d, point2_2d)
            #print("distance to line", distance_to_line)
            if distance_to_line <= distance_threshold:
                # Neuron is within the y-range and close enough to the plane
                number_neurons_on_implane += 1
    
    print("number of neurons on imaging plane", number_neurons_on_implane)
    return number_neurons_on_implane

#neurons_on_imaging_plane =neurons_on_imaging_plane(exp=4, pattern=0, mouse=0, amplitude=10, point1=[-9, 300, 16], angle_degrees=60)

def neurons_in_layer(exp, pattern, mouse, amplitude, point1= [-9, 300, 16], point2=[-9,300,189]):
    node_pos, spikes = get_spikes(exp, pattern, mouse, amplitude)
    #print("shape node pos", node_pos.shape)
    pos_filtered, spikes_filtered, threshold = filter_spikes(node_pos, spikes)
    #print("shape of pos_filtered is ", pos_filtered.shape)
    #print("spikes filtered", spikes_filtered)

    
    #Convert points to 2D (z, x) plane
    point1_2d = np.array([point1[2],point1[0]])
    point2_2d = np.array([point2[2],point2[0]])
    print("point1_2d", point1_2d)
    print("point2 2d", point2_2d)
    #print("shape point 1 2d is", point1_2d.shape)


    y_min = 100
    y_max = 430
    number_neurons_on_plane=0

    distance_threshold = 15.0

    for neuron in range(len(pos_filtered)):
        # Get x, y, z coordinates of the neuron
        x_neuron = pos_filtered[neuron, 0]
        y_neuron = pos_filtered[neuron, 1]
        z_neuron = pos_filtered[neuron, 2]

        # Check if y is within the range
        if y_min <= y_neuron <= y_max:
            neuron_2d = np.array([z_neuron, x_neuron])
            distance_to_line = point_line_distance(neuron_2d, point1_2d, point2_2d)
            #print("distance to line", distance_to_line)
            if distance_to_line <= distance_threshold:
                # Neuron is within the y-range and close enough to the plane
                number_neurons_on_plane += 1
    
    print("number of neurons on plane", number_neurons_on_plane)
    return number_neurons_on_plane

#neurons_layer =neurons_in_layer(exp=4, pattern=0, mouse=0, amplitude=10, point1=[-9, 300, 16],point2=[-9,300,189])

def plot_depth():
    plt.figure()
    markersize=7

    for mouse_i in [0,1,2]:
        print(mouse_i)
        n_60_10=neurons_on_imaging_plane(exp=4, pattern=0, mouse=mouse_i, amplitude=10, point1=[-9, 300, 16], angle_degrees=60)
        n_0_10=neurons_in_layer(exp=4, pattern=0, mouse=mouse_i, amplitude=10, point1=[-9, 300, 16],point2=[-9,300,189])
        n_60_20=neurons_on_imaging_plane(exp=4, pattern=0, mouse=mouse_i, amplitude=20, point1=[-9, 300, 16], angle_degrees=60)
        n_0_20=neurons_in_layer(exp=4, pattern=0, mouse=mouse_i, amplitude=20, point1=[-9, 300, 16],point2=[-9,300,189])


        
        plt.plot(1, n_0_10, 'o', markersize=markersize, color='orange')
        plt.plot(2, n_60_10, 'o', markersize=markersize, color='orange')
        plt.plot(1, n_0_20, 'o', markersize=markersize, color='blue')
        plt.plot(2, n_60_20, 'o', markersize=markersize, color='blue')

        # Connect dots with black lines
        plt.plot([1, 2], [n_0_10, n_60_10], color='black')
        plt.plot([1, 2], [n_0_20, n_60_20], color='black')

    # Plot "invisible" points to create a single legend entry for each color
    plt.plot([], [], 'o', markersize=markersize, color='blue', label='20 µA')
    plt.plot([], [], 'o', markersize=markersize, color='orange', label='10 µA')    

    plt.yticks(fontsize=15)
    plt.ylabel('Nb. of activated neurons', fontsize=20)
    plt.xticks([1,2], ['Imaging plane along layer', 'Imaging plane 60 degrees'], fontsize=15)
    #plt.yticks(np.arange(0, 21, 2), fontsize=15)
    #plt.ylim(0,20)
    plt.legend(fontsize=15)

    plt.show()

#plot_depth()





