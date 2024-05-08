import numpy as np
import matplotlib.pyplot as plt
import tifffile
import pickle
import os
import math
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import wilcoxon
from sklearn.neighbors import KernelDensity

class ExtractFluorescence():
    def __init__(self, directory, recording_name, nb_patterns):
        self.cell_data = None
        self.pattern_centroids = []
        self.nb_patterns = nb_patterns
        self.load_neurons(directory=directory, filename = recording_name + "_stimall_max_seg.npy")

        temp_fn = "temp_file_" + recording_name # Added extra "_" for non-downsampled file
        if os.path.exists(temp_fn):
            print("Loading in existing data for recording " + recording_name)
            self.read_intensities(temp_fn)
        else:
            print("Loading in new data for recording " + recording_name)
            self.load_recording(directory=directory, filename = recording_name + "_dff_not_downsampled_averaged.tif")
            self.get_intensities()
            self.normalize_intensities()
            self.write_intensities(temp_fn)

        self.get_cell_centroids()

    def load_neurons(self, directory, filename):
        self.neuron_fp = directory + "/" + filename
        try:
            array = np.load(self.neuron_fp, allow_pickle=True).item()
            self.masks = array["masks"]
        except FileNotFoundError:
            print("File not found at the specified path.")

        # Extract unique cell values from masks
        self.unique_cells = set(cell for row in self.masks for cell in row)
        self.nb_cells = len(self.unique_cells) - 1 # "0" is also seen as a cell but this is just the background without any cells

    def load_recording(self, directory, filename):
        self.recording_fp = directory + "/" + filename
        with tifffile.TiffFile(self.recording_fp) as tif:
            # Assuming the image stack is the first series in the TIFF file
            self.recording = tif.series[0].asarray()
            self.recording_length = self.recording.shape[0] # shape = (time, height, width)
        self.nb_pattern_timepoints = round(self.recording_length/self.nb_patterns) # 12 for a downsampled recording, 72 for a non-downsampled recording

    def write_intensities(self, fn):
        with open(fn, 'wb') as f:
            pickle.dump(self.cell_data, f)

    def read_intensities(self, fn):
        with open(fn, "rb") as f:
            self.cell_data = pickle.load(f)
        self.recording_length = self.cell_data[1]['intensities'].shape[0]
        self.nb_pattern_timepoints = round(self.recording_length/self.nb_patterns) # 12 for a downsampled recording, 72 for a non-downsampled recording

    def write_fluorescence_traces_denoising(self):
        # dF_traces.npy-files (neurons x time)
        # Best to include some mock data before the first peak, otherwise this is removed by the CASCADE denoising
        dF_traces = []
        for cell in self.cell_data.keys():
            dF_traces.append(self.cell_data[cell]['intensities'])
        dF_traces = np.array(dF_traces)
        # print(np.shape(dF_traces))
        np.save('dF_traces.npy', dF_traces)

    def read_fluorescence_traces_denoising(self):
        # Denoised with CASCADE (model = 'Global_EXC_10Hz_smoothing200ms_causalkernel')
        dF_traces = np.load('predictions_dF_traces.npy')
        print(np.shape(dF_traces))
        cells = [1, 2, 3, 4, 5]
        fig, axs = plt.subplots(len(cells), sharex=True)
        fig.subplots_adjust(hspace=0)
        ax_counter = 0
        for cell in cells:
            axs[ax_counter].plot(dF_traces[cell-1])
            ax_counter += 1
        fig.suptitle("Normalized fluorescence of different cells")
        plt.show()

    def get_intensities(self):
        # Initialize a dictionary to hold lists for each cell
        self.cell_data = {cell_value: {'intensities': []} for cell_value in self.unique_cells}

        # Loop over each time point
        for time_point in range(self.recording_length):
            # Initialize or reset sum of intensities and counts for this time point
            cell_sums = {cell_value: 0 for cell_value in self.unique_cells}
            cell_counts = {cell_value: 0 for cell_value in self.unique_cells}

            # Loop over each pixel and update sums and counts
            for y, row in enumerate(self.masks):
                for x, cell_value in enumerate(row):
                    if cell_value > 0: # Skip cell "0", this is just the background
                        intensity = self.recording[time_point, y, x]
                        cell_sums[cell_value] += intensity
                        cell_counts[cell_value] += 1

            # Compute and store the average for each cell at this time point
            for cell_value in self.unique_cells:
                if cell_counts[cell_value] > 0:
                    avg_intensity = cell_sums[cell_value] / cell_counts[cell_value]
                    self.cell_data[cell_value]['intensities'].append(avg_intensity)
                else:
                    self.cell_data[cell_value]['intensities'].append(None)
                    if cell_value > 0:
                        print("No intensities found for cell " + str(cell_value))

        self.cell_data.pop(0) # remove the background values from the dictionary

    def normalize_intensities(self):
        baseline_timepoints = [i for i in range(0, self.nb_pattern_timepoints*self.nb_patterns) if i % self.nb_pattern_timepoints > self.nb_pattern_timepoints/2] # Assume that the first 7 timepoints of each 12 are stimulation, and the last 5 of each 12 are baseline
        # print(baseline_timepoints)
        if len(baseline_timepoints) < 10:
            print("Very small (<10) number of timepoints for standard deviation computation. Please choose larger baseline_timepoints.")
        for cell in self.cell_data.keys():
            selected_intensities = [self.cell_data[cell]['intensities'][i] for i in baseline_timepoints if i < len(self.cell_data[cell]['intensities'])]
            self.cell_data[cell]['intensities'] -= np.mean(selected_intensities) # Subtract mean such that baseline (= pre-stim) fluorescence is zero
            self.cell_data[cell]['intensities'] /= max(self.cell_data[cell]['intensities']) # Normalize each cell between 0 and 1
            self.cell_data[cell]['std'] = np.std(selected_intensities) # Add standard deviation to the cell data

    def get_cell_centroids(self):
        # All centroids are in pixels, not in um
        if self.cell_data == None:
            print("Computing cell centroids without cell intensities.")
            self.cell_data = {cell_value: {'centroid': []} for cell_value in self.unique_cells}
            self.cell_data.pop(0)
        else:
            for cell in self.cell_data.keys():
                self.cell_data[cell]['centroid'] = []

        cell_coordinates = {}
        for y, row in enumerate(self.masks):
            for x, cell_value in enumerate(row):
                if cell_value not in cell_coordinates:
                    cell_coordinates[cell_value] = {"x_sum": 0, "y_sum": 0, "count": 0}

                # Accumulate the x and y coordinates
                cell_coordinates[cell_value]["x_sum"] += x
                cell_coordinates[cell_value]["y_sum"] += y
                cell_coordinates[cell_value]["count"] += 1

        # Calculate the central coordinates for each cell
        for cell_value, data in cell_coordinates.items():
            if cell_value > 0: # Ignore the background
                x_mean = data["x_sum"] / data["count"]
                y_mean = data["y_sum"] / data["count"]
                self.cell_data[cell_value]['centroid'] = [round(x_mean), round(y_mean)]

    def plot_centroids(self):
        x_values = [self.cell_data[cell]['centroid'][0] for cell in self.cell_data.keys()]
        y_values = [self.cell_data[cell]['centroid'][1] for cell in self.cell_data.keys()]

        plt.scatter(x_values, y_values, s=90)
        if len(self.pattern_centroids) > 0:
            for centroid in self.pattern_centroids:
                plt.scatter(centroid[0], centroid[1], color='red', marker='x', label='Centroid')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        plt.legend()
        plt.title('Centroids of cells')
        plt.show()

    def plot_fluorescence(self, cells = [1, 2, 3, 4, 5]):
        fig, axs = plt.subplots(len(cells), sharex=True)
        fig.subplots_adjust(hspace=0)
        ax_counter = 0
        for cell in cells:
            axs[ax_counter].plot(self.cell_data[cell]['intensities'])
            axs[ax_counter].hlines(3*self.cell_data[cell]['std'], 0, self.nb_pattern_timepoints*self.nb_patterns, linestyles='dashed')
            for pattern in range(0,self.nb_patterns):
                axs[ax_counter].plot(pattern*self.nb_pattern_timepoints, self.cell_data[cell]['intensities'][pattern*self.nb_pattern_timepoints], 'ro', markersize=3)
            ax_counter += 1
        fig.suptitle("Normalized fluorescence of different cells")
        plt.show()

    def get_weighted_centroids(self, patterns=None):
        # Returns: the (x, y) locations of the weighted centroids for each specified timepoint.
        if patterns == None:
            timepoints = [i*12 for i in range(0,self.nb_patterns)]
        elif max(patterns) <= self.nb_patterns:
            timepoints = [i*12 for i in [pattern-1 for pattern in patterns]] # start counting patterns from 1
        else:
            print('Too large pattern given, max pattern is ' + str(self.nb_patterns))
        for timepoint in timepoints:
            total_weight = 0
            weighted_sum_x = 0
            weighted_sum_y = 0
            for cell in self.cell_data.keys():
                cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
                if cell_intensity > 3*self.cell_data[cell]['std']:
                    cell_centroid = self.cell_data[cell]['centroid']
                    weighted_sum_x += cell_centroid[0] * cell_intensity
                    weighted_sum_y += cell_centroid[1] * cell_intensity
                    total_weight += cell_intensity
            if total_weight > 0:
                weighted_centroid_x = weighted_sum_x / total_weight
                weighted_centroid_y = weighted_sum_y / total_weight

            self.pattern_centroids.append([weighted_centroid_x, weighted_centroid_y])

class ProcessFluorescence():
    def __init__(self, cell_data, nb_patterns):
        self.cell_data = cell_data
        self.nb_patterns = nb_patterns
        self.recording_length = self.cell_data[1]['intensities'].shape[0]
        self.nb_pattern_timepoints = round(self.recording_length/self.nb_patterns) # 12 for a downsampled recording, 72 for a non-downsampled recording
        self.get_common_intensities()

    def plot_centroids(self):
        pass

    def get_common_intensities(self, patterns=None):
        self.common_intensities = {cell_value: 0 for cell_value in self.cell_data.keys()}
        if patterns == None: # Look at all patterns
            timepoints = [i*12 for i in range(0,self.nb_patterns)]
        elif max(patterns) <= self.nb_patterns:
            timepoints = [i*12 for i in [pattern-1 for pattern in patterns]] # start counting patterns from 1
        else:
            print('Too large pattern given, max pattern is ' + str(self.nb_patterns))
        for cell in self.cell_data.keys():
            active_intensities = [max(0, self.cell_data[cell]['intensities'][timepoint]) for timepoint in timepoints] # xxxxxxxxxx # if self.cell_data[cell]['intensities'][timepoint] > 3*self.cell_data[cell]['std']
            self.common_intensities[cell] = min(active_intensities)
        # print(self.common_intensities)

    def get_weighted_centroid(self, pattern=1):
        # Centroid weighted by fluorescence
        # Returns: the (x, y) location of the weighted centroids for a specified pattern
        timepoint = self.nb_pattern_timepoints * (pattern-1)
        total_weight = 0
        weighted_sum_x = 0
        weighted_sum_y = 0
        unique_total_weight = 0
        unique_weighted_sum_x = 0
        unique_weighted_sum_y = 0
        for cell in self.cell_data.keys():
            cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
            if cell_intensity > 3*self.cell_data[cell]['std']:
                unique_cell_intensity = cell_intensity - self.common_intensities[cell]
                cell_centroid = self.cell_data[cell]['centroid']
                weighted_sum_x += cell_centroid[0] * cell_intensity
                weighted_sum_y += cell_centroid[1] * cell_intensity
                total_weight += cell_intensity

                unique_weighted_sum_x += cell_centroid[0] * unique_cell_intensity
                unique_weighted_sum_y += cell_centroid[1] * unique_cell_intensity
                unique_total_weight += unique_cell_intensity

        if total_weight > 0:
            weighted_centroid_x = weighted_sum_x / total_weight
            weighted_centroid_y = weighted_sum_y / total_weight
            unique_weighted_centroid_x = unique_weighted_sum_x / unique_total_weight
            unique_weighted_centroid_y = unique_weighted_sum_y / unique_total_weight
            return [weighted_centroid_x, weighted_centroid_y], [unique_weighted_centroid_x, unique_weighted_centroid_y]
        else:
            return None

    def fluorescence_direction(self, pattern=1, return_elec = [0,0], neurons='all'):
        central_elec = [241, 155]                        # [x, y] of electrode, [241, 155] for Slice4 from 31/01
        centroid, unique_centroid = self.get_weighted_centroid(pattern=pattern)
        # current_direction =
        # neural_direction =
        self.plot_active_neurons(pattern=pattern, neurons=neurons)
        plt.scatter(centroid[0], centroid[1], color='red', marker='x', label='Centroid')
        plt.scatter(unique_centroid[0], unique_centroid[1], color='purple', marker='x', label='Unique centroid')
        plt.plot([central_elec[0], return_elec[0]],[central_elec[1], return_elec[1]], 'green', linestyle="--")
        plt.scatter(central_elec[0], central_elec[1], color='green', marker='x', label='Central elec')
        plt.scatter(return_elec[0], return_elec[1], color='green', marker='x', label='Return elec')

        plt.legend()
        plt.title('Active neurons (' + neurons + ') for stim. pattern ' + str(pattern))
        # plt.show()

    def discriminate_signed_rank(self, pattern1=1, pattern2=None):
        '''
        Use the Wilcoxon signed-rank test to get a p-value as index of separability between the two neuronal populations.
        '''
        timepoint1 = self.nb_pattern_timepoints * (pattern1-1)
        if pattern2 != None:
            timepoint2 = self.nb_pattern_timepoints * (pattern2-1)
        else:
            timepoint2 = None

        fluorescence1 = []
        fluorescence2 = []
        color_value = []

        plotting_gradient = 'X'
        for cell in self.cell_data.keys():
            cell_intensity1 = self.cell_data[cell]['intensities'][timepoint1] # intensity at stimulation point 1
            cell_intensity2 = self.cell_data[cell]['intensities'][timepoint2] if timepoint2 != None else timepoint1 # intensity at stimulation point 2
            if cell_intensity1 > 3*self.cell_data[cell]['std'] or cell_intensity2 > 3*self.cell_data[cell]['std']:
                if pattern2 == None:
                    cell_intensity2 = self.common_intensities[cell]
                # todo: use fluorescence of non-active cell or use 0?
                cell_centroid = self.cell_data[cell]['centroid']
                fluorescence1.append(cell_intensity1)
                fluorescence2.append(cell_intensity2)
                if plotting_gradient == 'X':
                    color_value.append(cell_centroid[0]/397) # Look at X-gradient
                elif plotting_gradient == 'Y':
                    color_value.append(cell_centroid[1]/380) # Look at X-gradient

        color_value='blue' # If no gradient should be plotted. Comment out if gradient should be plotted

        plt.figure()
        plt.scatter(fluorescence1, fluorescence2, c=color_value, s=30, cmap='viridis')
        plt.xlabel('Fluorescence response for pattern ' + str(pattern1))
        plt.ylabel('Fluorescence response for pattern ' + str(pattern2))
        plt.title('Correlation of fluorescence values')

        signed_rank = wilcoxon(fluorescence1, fluorescence2) # Apply Wilcoxon test
        print('P-value for Wilcoxon signed-rank test for stim patterns ' + str(pattern1) + ' and ' + str(pattern2) + ' is ' + str(round(signed_rank.pvalue,5)))
        return(signed_rank.pvalue)

    def discriminate_svm(self, pattern1=1, pattern2=None):
        '''
        Build an SVM
        Usually (advice from Dylan): n-dimensional space with m samples. n = number of neurons (e.g. 80), m = number of repetitions (e.g. 20).
        In this case: only 1 to 6 samples in an 80-dimensional space, would be "too easy" to separate.
        Better: 3-dimensional space (X, Y, fluorescence) but then SVM doesn't really make sense anymore
        '''
        timepoint1 = self.nb_pattern_timepoints * (pattern1-1)
        if pattern2 != None:
            timepoint2 = self.nb_pattern_timepoints * (pattern2-1)
        else:
            timepoint2 = None

        data = []
        labels = []
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for cell in self.cell_data.keys():
            cell_intensity1 = self.cell_data[cell]['intensities'][timepoint1] # intensity at stimulation point 1
            cell_intensity2 = self.cell_data[cell]['intensities'][timepoint2] if timepoint2 != None else timepoint1 # intensity at stimulation point 2
            if cell_intensity1 > 3*self.cell_data[cell]['std'] or cell_intensity2 > 3*self.cell_data[cell]['std']:
                if pattern2 == None:
                    cell_intensity2 = self.common_intensities[cell] # If pattern2 == None, compare the unique neurons with the common neurons
                cell_centroid = self.cell_data[cell]['centroid']
                # For each cell, add 2 samples: 1 for pattern 1, and 1 for pattern 2
                data.append([cell_centroid[0], cell_centroid[1], cell_intensity1])
                labels.append(0)
                data.append([cell_centroid[0], cell_centroid[1], cell_intensity2])
                labels.append(1)
                ax.scatter(cell_centroid[0], cell_centroid[1], cell_intensity1, c='blue')
                ax.scatter(cell_centroid[0], cell_centroid[1], cell_intensity2, c='red')

        data = np.array(data)
        labels = np.array(labels)

        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_zlabel('Fluorescence')

        data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.4, random_state=42) # Divide in training and test set
        clf = svm.SVC(kernel="rbf", C=1, random_state=42) # Choose model for SVM
        clf = make_pipeline(StandardScaler(), clf) # Add normalization automatically to the beginning of the SVM
        clf.fit(data_train, labels_train) # Train the SVM
        score = clf.score(data_test, labels_test) # Test the SVM
        return score

    def get_overlap(self, pattern1=1, pattern2=None):
        # Get number of neurons for a certain pattern, and overlap between neurons for different patterns
        # If pattern2 == None, compare the unique neurons with the common neurons

        timepoint1 = self.nb_pattern_timepoints * (pattern1-1)
        if pattern2 != None:
            timepoint2 = self.nb_pattern_timepoints * (pattern2-1)
        else:
            timepoint2 = None

        x = []
        y = []
        overlap = []
        alpha = []
        for cell in self.cell_data.keys():
            cell_intensity1 = self.cell_data[cell]['intensities'][timepoint1] # intensity at stimulation point 1
            cell_intensity2 = self.cell_data[cell]['intensities'][timepoint2] if timepoint2 != None else timepoint1 # intensity at stimulation point 2
            if cell_intensity1 > 3*self.cell_data[cell]['std'] or cell_intensity2 > 3*self.cell_data[cell]['std']:
                if pattern2 == None:
                    cell_intensity2 = self.common_intensities[cell]

                cell_centroid = self.cell_data[cell]['centroid']
                # print(cell_centroid, cell_intensity1, cell_intensity2)
                x.append(cell_centroid[0])
                y.append(cell_centroid[1])
                overlap.append(cell_intensity1-cell_intensity2)
                alpha.append(max(cell_intensity1,cell_intensity2))
        plt.figure()
        plt.scatter(x, y, s=90, c=overlap, alpha=alpha, cmap='rainbow')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        cbar = plt.colorbar(ticks=[min(overlap), 0, max(overlap)])
        if pattern2 == None:
            pattern2 = 'Common intensity'
        cbar.ax.set_yticklabels(['Only by ' + str(pattern2), 'Equal', 'Only by ' + str(pattern1)])
        # plt.legend()
        plt.title('Overlapping activity for stim. patterns ' + str(pattern1) + ' & ' + str(pattern2))
        # plt.show()

        # Calculate absolute differences for each coordinate pair
        absolute_overlap = np.abs(overlap)
        # Normalize the absolute differences to range [0, 1]
        normalized_overlap = absolute_overlap / np.max(absolute_overlap)
        # Compute the mean of the normalized absolute differences
        similarity_metric = np.mean(normalized_overlap)

        fig = plt.figure()
        plt.hist(normalized_overlap)
        xt = plt.xticks()
        loc0 = np.where(xt[0]==0)[0][0]
        loc1 = np.where(xt[0]>0.99)[0][0] # For some reason ==1 is not recognized
        label0 = xt[1][loc0]
        label0.set_text(label0.get_text()+'\nSimilar')
        label1 = xt[1][loc1]
        label1.set_text(label1.get_text()+'\nDifferent')
        xt[1][loc0] = label0
        xt[1][loc1] = label1
        plt.xticks(xt[0][loc0:loc1+1], xt[1][loc0:loc1+1])
        plt.xlabel('Difference in fluorescence')
        plt.ylabel('Weighted cell count')
        plt.title('Similarity of pattern ' + str(pattern1) + ' & ' + str(pattern2) + ': ' + str(np.round(similarity_metric,2)))
        # plt.show()

        return x, y, overlap, similarity_metric

    def kernel_density_estimate(self, pattern=1):
        '''
        2D Kernel Density Estimate of the data
        '''
        coordinates = []
        fluorescence = []
        timepoint = self.nb_pattern_timepoints * (pattern-1)
        for cell in self.cell_data.keys():
            cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
            if cell_intensity > 3*self.cell_data[cell]['std']:
                # cell_intensity -= self.common_intensities[cell]
                coordinates.append(self.cell_data[cell]['centroid'])
                fluorescence.append(cell_intensity)
        coordinates = np.array(coordinates)

        kde = KernelDensity(bandwidth=100, kernel='gaussian') # Choose model and parameters
        kde.fit(coordinates, sample_weight=fluorescence) # Train model

        grid_size = 100 # 100 points in x and in y direction
        x_grid, y_grid = np.meshgrid(np.linspace(0, 397, grid_size), np.linspace(0, 380, grid_size))
        grid_points = np.vstack([x_grid.ravel(), y_grid.ravel()]).T
        density = np.exp(kde.score_samples(grid_points)).reshape(x_grid.shape) # Evaluate model for all points on the grid

        fig = plt.figure()
        plt.pcolormesh(x_grid, y_grid, density, shading='auto')
        plt.scatter(coordinates[:,0], coordinates[:,1], c=fluorescence, cmap='viridis', edgecolors='k', linewidths=1)
        plt.colorbar(label='Values')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        # plt.legend()
        plt.title('Kernel Density Estimate for stim. pattern ' + str(pattern))
        # plt.show()
        plt.close()

        return coordinates, fluorescence, x_grid, y_grid, density

    def projected_kernel_density_estimate(self, pattern=1, point1=(0,0), point2=(397,380)):
        '''
        Projection of the data before the Kernel Density Estimate is useful when trying the understand the spatial
        distribution of cell activity along a specific axis. The density estimate now reflects the distribution
        of the data along the projected data.
        If projection would only take place after the Kernel Density Estimate, then you just integrate the density
        function, but the density estimate is still based on a 2D-distribution.
        '''
        direction_vector = np.array([point2[0] - point1[0], point2[1] - point1[1]])
        projected_points = []
        fluorescence = []
        timepoint = self.nb_pattern_timepoints * (pattern-1)

        fig, ax = plt.subplots(3, 1, figsize=(5, 15))
        ax = ax.flatten()

        for cell in self.cell_data.keys():
            cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
            if cell_intensity > 3*self.cell_data[cell]['std']:
                # cell_intensity -= self.common_intensities[cell]
                centroid = self.cell_data[cell]['centroid']

                # Project the cell on the user-defined axis
                t = np.dot(np.array(centroid) - np.array(point1), direction_vector) / np.dot(direction_vector, direction_vector)
                projected_points.append(t)
                projected_coordinate = [point1[0] + t * direction_vector[0], point1[1] + t * direction_vector[1]]
                fluorescence.append(cell_intensity)
                ax[0].scatter(centroid[0], centroid[1], s=90, c="blue", alpha=cell_intensity)
                ax[0].scatter(projected_coordinate[0], projected_coordinate[1], s=10, c="red", alpha=cell_intensity)
        projected_points = np.array(projected_points)
        fluorescence = np.array(fluorescence)

        ax[0].set_xlabel('X Coordinate')
        ax[0].set_ylabel('Y Coordinate')
        ax[0].set_xlim([0, 397])
        ax[0].set_ylim([0, 380])
        ax[0].invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        ax[0].set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal

        # Plot for original points, line, and projected points
        ax[0].plot([point1[0], point2[0]], [point1[1], point2[1]], color='blue', label='Line')
        ax[0].set_title('Projection of Points onto Line')
        ax[0].legend()

        # Plot for points along the projection axis
        ax[1].scatter(projected_points, fluorescence, color='blue')
        ax[1].set_xticks([])
        ax[1].set_ylabel('Fluorescence')
        ax[1].set_title('Points along projection Axis')
        ax[1].set_aspect('auto')

        # Define grid along the projected axis
        grid_size = 100
        grid = np.linspace(min(projected_points), max(projected_points), grid_size).reshape(-1, 1)

        # Perform kernel density estimation
        kde = KernelDensity(bandwidth=0.3, kernel='gaussian') # Grid only goes from about -2 to +2 so take smaller bandwidth than 2D-KDE where grid goes from 0 to 397

        kde.fit(projected_points.reshape(-1, 1), sample_weight=fluorescence)
        density = np.exp(kde.score_samples(grid))

        # Plot 1D density function along the projected axis on ax1
        ax[2].plot(grid, density, color='red', linestyle='-')
        ax[2].set_xlabel('Distance along the projected axis')
        ax[2].set_ylabel('Density')
        ax[2].set_title('1D Kernel Density Estimate along Projected Axis')

        ax[2].set_xticks([0,1],["Point 1", "Point 2"])
        ax[2].sharex(ax[1])

        plt.tight_layout()
        # plt.show()
        plt.close()

        return grid, density

    def full_kde(self, pattern=1):
        # For slice 4 - 0:
        # central elec: [241, 155]
        # 1) elec 1: [157, 60]
        # 2) elec 9: [328, 254]
        # 3) elec 3,4: [12, 75]
        # 4) elec 5: [142, 224]
        # 5) elec 8: [415, 353]

        central_elec = (241, 155)
        return_elec_layer = (328, 254)
        return_elec_col = (142, 224)
        return_elecs = [(157, 60), (328, 254), (12, 75), (142, 224), (415, 353)]
        return_elec = return_elecs[pattern-1]

        coordinates, fluorescence, x_grid, y_grid, xy_density = self.kernel_density_estimate(pattern=pattern)
        proj_layer_grid, proj_layer_density = self.projected_kernel_density_estimate(pattern=pattern, point1=central_elec, point2=return_elec_layer)
        max_layer = proj_layer_grid[np.argmax(proj_layer_density)][0] # Compute the grid value that corresponds to the index of the highest density value
        max_layer_xy = (central_elec[0] + max_layer*(return_elec_layer[0]-central_elec[0]), central_elec[1] + max_layer*(return_elec_layer[1]-central_elec[1])) # Compute the X and Y for this maximal value
        proj_col_grid, proj_col_density = self.projected_kernel_density_estimate(pattern=pattern, point1=central_elec, point2=return_elec_col)
        max_col = proj_col_grid[np.argmax(proj_col_density)][0]
        max_col_xy = (central_elec[0] + max_col*(return_elec_col[0]-central_elec[0]), central_elec[1] + max_col*(return_elec_col[1]-central_elec[1]))

        # plt.close('all') # Close all other figures

        fig = plt.figure(figsize=(8,12))

        ax1 = plt.subplot2grid((4, 2), (0, 0), colspan=1, rowspan=2)  # 1st row, 1st column, spanning 1 column
        ax2 = plt.subplot2grid((4, 2), (0, 1), colspan=1, rowspan=2)  # 1st row, 2nd column, spanning 1 column
        ax3 = plt.subplot2grid((4, 2), (2, 0), colspan=2)  # 2nd row, 1st column, spanning 2 columns
        ax4 = plt.subplot2grid((4, 2), (3, 0), colspan=2)  # 3rd row, 1st column, spanning 2 columns

        # ax1.plot([central_elec[0], return_elec_layer[0]],[central_elec[1], return_elec_layer[1]], 'limegreen', linestyle="--", label='Along layer IV')
        ax1.axline(central_elec, return_elec_layer, color='limegreen', label='Along layer')
        # ax1.plot([central_elec[0], return_elec_col[0]],[central_elec[1], return_elec_col[1]], 'darkgreen', linestyle="--", label='Along cortical column')
        ax1.axline(central_elec, return_elec_col, color='darkgreen', label='Along column')
        ax1.scatter(coordinates[:,0], coordinates[:,1], s=90, c="blue", alpha=fluorescence)
        ax1.scatter(central_elec[0], central_elec[1], color='orange', s=110, marker='s', label='Central elec. (L4)', zorder=3)
        # ax1.scatter(return_elec_col[0], return_elec_col[1], color='darkgreen', label='Adjacent elec. (L2/3)')
        # ax1.scatter(return_elec_layer[0], return_elec_layer[1], color='limegreen', label='Adjacent elec. (L4)')
        ax1.scatter(return_elec[0], return_elec[1], color='gold', s=110, marker='s', label='Return elec.', zorder=3)
        ax1.scatter(max_layer_xy[0], max_layer_xy[1], color='red', marker='*', s=120, label='Max density', zorder=3)
        ax1.scatter(max_col_xy[0], max_col_xy[1], color='red', marker='*', s=120, zorder=3)

        ax1.set_xlabel('X Coordinate')
        ax1.set_ylabel('Y Coordinate')
        ax1.set_xlim([0, 397])
        ax1.set_ylim([0, 380])
        ax1.invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        ax1.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        ax1.legend(fontsize='8', loc='upper right')

        pcm = ax2.pcolormesh(x_grid, y_grid, xy_density, shading='auto')
        ax2.scatter(coordinates[:,0], coordinates[:,1], c=fluorescence, cmap='viridis', edgecolors='k', linewidths=1)
        fig.colorbar(pcm, ax=ax2, label='Values')
        ax2.set_xlabel('X Coordinate')
        ax2.set_ylabel('Y Coordinate')
        ax2.set_xlim([0, 397])
        ax2.set_ylim([0, 380])
        ax2.invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        ax2.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        ax2.set_title('2D Kernel Density Estimate')

        ax3.plot(proj_layer_grid, proj_layer_density, color='red', linestyle='-')
        ax3.set_xlabel('Distance along layer IV')
        ax3.set_ylabel('Density')
        ax3.set_xticks([0,1],["Central elec. (L4)", "Adjacent elec. (L4)"])
        ax3.set_title('1D Kernel Density Estimates along Projected Axes')

        ax4.plot(proj_col_grid, proj_col_density, color='red', linestyle='-')
        ax4.set_xlabel('Distance along the cortical column')
        ax4.set_ylabel('Density')
        ax4.set_xticks([0,1],["Central elec. (L4)", "Adjacent elec. (L2/3)"])
        # ax4.set_title('1D Kernel Density Estimate along Projected Axis')

        fig.suptitle('Kernel Density Estimate for stimulation pattern ' + str(pattern))
        plt.tight_layout(h_pad=4)
        # plt.show()

        return max_layer_xy, max_col_xy

    def fluorescence_distance(self, pattern=1, return_elec = [0,0]):
        # fluorescence strength in function of direction of current/distance to return electrode
        bin_size = 15 # number of pixels per bin = radius of concentric circle around the electrode
        intensities = [0]*(round(400/bin_size))    # divide in bins of bin_size pixels around electrode - important! is in pixels, not in um! Entire FOV = +- 400
        central_elec = [241, 155]                        # [x, y] of electrode, [241, 155] for Slice4 from 31/01
        elec = [np.mean([central_elec[0], return_elec[0]]), np.mean([central_elec[1], return_elec[1]])]
        print("Calculating distance to " + str(elec))
        timepoint = self.nb_pattern_timepoints * (pattern-1)
        for cell in self.cell_data.keys():
            cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
            if cell_intensity > 3*self.cell_data[cell]['std']:
                # cell_intensity -= self.common_intensities[cell]
                cell_centroid = self.cell_data[cell]['centroid']
                distance = math.sqrt((cell_centroid[0]-elec[0])**2 + (cell_centroid[1]-elec[1])**2)
                print(np.ceil(distance/bin_size))
                intensities[round(np.ceil(distance/bin_size))] += cell_intensity/(np.ceil(distance/bin_size)**2)
        # plt.plot(intensities, label=str(pattern) + " to " + str(elec))
        plt.plot(self.moving_average(intensities[1:],5), label=str(pattern) + " to " + str(elec))
        plt.legend()
        plt.xlabel('Distance to virtual electrode (nb of 15um bins)')
        plt.ylabel('Weighted cell count')
        plt.title('Activity in function of the distance')
        # plt.show()

    def moving_average(self, arr, window_size):
        ma = []
        i = 0
        while i < len(arr) - window_size + 1:
            ma.append(np.mean(arr[i:i+window_size]))
            i += 1
        return ma

    def plot_common_neurons(self):
        plt.figure()
        for cell in self.cell_data.keys():
            cell_intensity = self.common_intensities[cell]
            cell_centroid = self.cell_data[cell]['centroid']
            plt.scatter(cell_centroid[0], cell_centroid[1], s=90, c="blue", alpha=cell_intensity)
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        plt.legend()
        plt.title('Common neurons')
        # plt.show()

    def plot_active_neurons(self, pattern=1, neurons='unique'):
        # Only do this for "unique subset of neurons", otherwise too many neurons are taken into account (too many neurons are slightly active, weighted version would also work)
        timepoint = self.nb_pattern_timepoints * (pattern-1)
        plt.figure()
        for cell in self.cell_data.keys():
            cell_intensity = self.cell_data[cell]['intensities'][timepoint] # intensity at stimulation point
            if cell_intensity > 3*self.cell_data[cell]['std']:
                cell_centroid = self.cell_data[cell]['centroid']
                if neurons == 'all':
                    cell_intensity = cell_intensity
                elif neurons == 'unique':
                    cell_intensity = cell_intensity - self.common_intensities[cell]
                plt.scatter(cell_centroid[0], cell_centroid[1], s=90, c="blue", alpha=cell_intensity)
                # print(str(cell_intensity) + ', ', end="")
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        plt.legend()
        plt.title('Active neurons (' + neurons + ') for stim. pattern ' + str(pattern))
        # plt.show()

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment0 = ExtractFluorescence(directory=directory, recording_name="Slice4 - 0", nb_patterns=5)
# experiment8 = ExtractFluorescence(directory=directory, recording_name="Slice4 - 8", nb_patterns=6)
comparison0 = ProcessFluorescence(experiment0.cell_data, nb_patterns=5)
# comparison8 = ProcessFluorescence(experiment8.cell_data, nb_patterns=6)

# Slice 4 - 0
# central elec: [241, 155]
# 1) elec 1: [157, 60]
# 2) elec 9: [328, 254]
# 3) elec 3,4: [12, 75]
# 4) elec 5: [142, 224]
# 5) elec 8: [415, 353]
return_elecs = [[157, 60], [328, 254], [12, 75], [142, 224], [415, 353]]

pattern=2
comparison0.fluorescence_direction(pattern=pattern, return_elec=return_elecs[pattern-1], neurons='unique')
#comparison0.kernel_density_estimate(pattern=pattern)
#comparison0.discriminate_svm(pattern1=2, pattern2=4)
#plt.show()

comparison0.discriminate_signed_rank(pattern1=1, pattern2=3)
plt.show()

exit()


#comparison0.kernel_density_estimate(pattern=4)
#comparison0.projected_kernel_density_estimate(pattern=4, point1=(100, 100), point2=(300, 300))

pattern=5
comparison0.full_kde(pattern=pattern)
plt.show()

exit()

for pattern in [2, 4]:
    comparison0.fluorescence_direction(pattern=pattern, return_elec=return_elecs[pattern-1], neurons='unique')
comparison0.get_overlap(pattern1=2, pattern2=4)
print(comparison0.discriminate_signed_rank(pattern1=1, pattern2=3))
comparison0.kernel_density_estimate(pattern=4)
comparison0.projected_kernel_density_estimate(pattern=2, point1=(1, 100), point2=(200, 300))




plt.show()
exit()

'''
comparison0.plot_common_neurons()
comparison0.plot_active_neurons(pattern=4,neurons='all')
comparison0.plot_active_neurons(pattern=4,neurons='unique')
comparison0.plot_active_neurons(pattern=2,neurons='all')
comparison0.plot_active_neurons(pattern=2,neurons='unique')
for pattern in range(1,6):
    comparison0.fluorescence_direction(pattern=pattern, return_elec=return_elecs[pattern-1], neurons='unique')

comparison0.get_overlap(pattern1=2, pattern2=4)
'''

'''
plt.figure()
for pattern in range(1,6):
    comparison0.fluorescence_distance(pattern=pattern, return_elec=return_elecs[pattern-1])
plt.show()
'''

'''
plt.figure()
comparison0.fluorescence_distance(pattern=2, return_elec = [142, 224])
comparison0.fluorescence_distance(pattern=2, return_elec = [328, 254])
comparison0.fluorescence_distance(pattern=4, return_elec = [142, 224])
comparison0.fluorescence_distance(pattern=4, return_elec = [328, 254])
plt.show()

exit()

plt.figure()
comparison8.fluorescence_distance(pattern=2, return_elec = [142, 224])
comparison8.fluorescence_distance(pattern=2, return_elec = [328, 254])
# comparison8.fluorescence_distance(pattern=1, return_elec = [142, 224])
# comparison8.fluorescence_distance(pattern=1, return_elec = [328, 254])
plt.show()
'''