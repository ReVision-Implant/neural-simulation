import numpy as np
import matplotlib.pyplot as plt
import tifffile
import pickle
import os
import math
from scipy.stats import wilcoxon, pearsonr, spearmanr
from scipy.optimize import curve_fit
from sklearn.neighbors import KernelDensity
from matplotlib.patches import Ellipse, FancyArrowPatch

class ExtractFluorescence():
    '''
    Structure of the cell_data from a recording: dictionary with all the cells: cell_data[1] for cell 1, etc
    Structure of the cell_data[x]: dictionary with all the information for this cell:
    'location': (x,y) of the center of the cell
    'intensities': time traces of the intensity for each pattern: [[pattern 1],[pattern 2], etc]. The first time sample is the sample right after stimulation
    'peaks': intensity peaks for each pattern: [pattern 1, pattern 2, etc]
    'active': mask that determines whether a cell is active for a certain pattern (1=active, 0=inactive): [pattern 1, pattern 2, etc]
    'std': standard deviation of the pre-stimulation
    '''
    def __init__(self, directory, recording_name, nb_patterns):
        self.cell_data = None
        self.nb_patterns = nb_patterns
        self.cells_to_pop = [0]

        processed_fn = "cell_data_" + recording_name #
        if os.path.exists(directory + "/" + processed_fn):
            print("Loading in existing data for recording " + recording_name)
            self.read_cell_data(directory, processed_fn)
        else:
            print("Loading in new data for recording " + recording_name)
            self.load_neurons(directory=directory, filename = recording_name + "_stimall_seg.npy") # Load in the masks for the neurons
            self.load_recording(directory=directory, filename = recording_name + "_dff.tif") # can also be with "_downsampled" etc behind it
            self.get_cell_intensities()
            self.split_cell_intensities()
            self.get_cell_peaks()
            self.get_cell_std()
            self.get_cell_centroids()
            self.write_cell_data(directory, processed_fn)
            print("Finished writing new data for recording " + recording_name)

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
        if self.nb_pattern_timepoints != 72:
            print('Possibly something wrong with the number of patterns. Should be : ' + str(self.recording_length/self.nb_pattern_timepoints))

    def get_cell_intensities(self):
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
                        if intensity >= 0: # Filter out the NaN (not clear why this sometimes happens)
                            cell_sums[cell_value] += intensity
                            cell_counts[cell_value] += 1
                        else:
                            # print(str(intensity) + ' intensity filtered out for cell ' + str(cell_value) + ' at x=' + str(x) + ' and y=' + str(y))
                            if cell_value not in self.cells_to_pop:
                                self.cells_to_pop.append(cell_value)
                                print('Added cell ' + str(cell_value) + ' to cells_to_pop because of NaN values.')

            # Compute and store the average for each cell at this time point
            for cell_value in self.unique_cells:
                if cell_counts[cell_value] > 0:
                    avg_intensity = cell_sums[cell_value] / cell_counts[cell_value]
                    self.cell_data[cell_value]['intensities'].append(avg_intensity)
                else:
                    self.cell_data[cell_value]['intensities'].append(None)
                    if cell_value > 0:
                        print("No intensities found for cell " + str(cell_value))

        for cell in self.cells_to_pop:
            self.cell_data.pop(cell) # remove the background values from the dictionary

    def split_cell_intensities(self):
        for cell in self.cell_data.keys():
            intensities = [[]]*self.nb_patterns
            for pattern in range(self.nb_patterns):
                intensities[pattern] = self.cell_data[cell]['intensities'][pattern*self.nb_pattern_timepoints:(pattern+1)*self.nb_pattern_timepoints]
            self.cell_data[cell]['intensities'] = intensities

    def get_cell_peaks(self):
        for cell in self.cell_data.keys():
            peaks = []
            peak_locations = []
            for pattern in range(self.nb_patterns):
                pattern_peak = max(self.cell_data[cell]['intensities'][pattern][:int(self.nb_pattern_timepoints/12*1.6)]) # Look for the maximum in the first second after stimulation (equivalent to 1.6 samples for downsampled data)
                pattern_peak_loc = self.cell_data[cell]['intensities'][pattern].index(pattern_peak)
                peaks.append(pattern_peak)
                peak_locations.append(pattern_peak_loc)
            self.cell_data[cell]['peaks'] = peaks
            self.cell_data[cell]['peak_locations'] = peak_locations

    def get_cell_std(self):
        for cell in self.cell_data.keys():
            background = []
            for pattern in range(self.nb_patterns):
                background.extend(self.cell_data[cell]['intensities'][pattern][-int(self.nb_pattern_timepoints/12*5):]) # Assume that the first 7 timepoints of each 12 are stimulation, and the last 5 of each 12 are baseline
            mean = np.mean(background)
            std = np.std(background)
            active = []
            for pattern in range(self.nb_patterns):
                if self.cell_data[cell]['peaks'][pattern] > mean + 3*std:
                    active.append(1)
                else:
                    active.append(0)
            self.cell_data[cell]['mean'] = mean
            self.cell_data[cell]['std'] = std
            self.cell_data[cell]['active'] = active

    def get_cell_centroids(self):
        # All centroids are in pixels, not in um todo change to um
        if self.cell_data == None:
            print("Warning: computing cell centroids without cell intensities.")
            self.cell_data = {cell_value: {'centroid': []} for cell_value in self.unique_cells}
            for cell in self.cells_to_pop:
                self.cell_data.pop(cell) # remove the background values from the dictionary
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
            if cell_value not in self.cells_to_pop: # Ignore the background
                x_mean = data["x_sum"] / data["count"]
                y_mean = data["y_sum"] / data["count"]
                self.cell_data[cell_value]['centroid'] = [round(x_mean), round(y_mean)]

    def write_cell_data(self, directory, fn):
        with open(directory + "/" + fn, 'wb') as f:
            pickle.dump(self.cell_data, f)

    def read_cell_data(self, directory, fn):
        with open(directory + "/" + fn, "rb") as f:
            self.cell_data = pickle.load(f)
        self.recording_length = len(self.cell_data[next(iter(self.cell_data.keys()))]['intensities'][0])
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


class ProcessFluorescence():
    def __init__(self, cell_data, electrode_locations):
        self.cell_data = cell_data
        self.nb_patterns = len(self.cell_data[next(iter(self.cell_data.keys()))]['intensities'])
        self.nb_pattern_timepoints = len(self.cell_data[next(iter(self.cell_data.keys()))]['intensities'][0])
        self.get_common_activity()
        self.max_activity = np.max([self.cell_data[cell]['peaks'] for cell in self.cell_data.keys()]) # Only for visualization purposes, alpha has to be between 0 and 1
        self.load_electrode_locations(electrode_locations[0], electrode_locations[1], electrode_locations[2], electrode_locations[3])

    def load_electrode_locations(self, central_electrode, return_electrodes, return_electrode_column, return_electrode_layer):
        # Electrode coordinates
        # For slice 4 - 0 (31/01/2024):
        # central elec: [241, 155]
        # 0) elec 1: [157, 60]
        # 1) elec 9: [328, 254]
        # 2) elec 3,4: [12, 75]
        # 3) elec 5: [142, 224]
        # 4) elec 8: [415, 353]
        # self.central_electrode = (241, 155)
        # self.return_electrodes = [(157, 60), (328, 254), (12, 75), (142, 224), (415, 353)]
        # self.return_electrode_layer = (328, 254)
        # self.return_electrode_column = (142, 224)
        if not len(return_electrodes) == self.nb_patterns:
            print('Give the correct number of return electrodes: ' + str(self.nb_patterns))
            return
        self.central_electrode = central_electrode
        self.return_electrodes = return_electrodes
        self.return_electrode_layer = return_electrode_layer
        self.return_electrode_column = return_electrode_column

    def get_common_activity(self, patterns=None):
        for cell in self.cell_data.keys():
            self.cell_data[cell]['common_activity'] = min([self.cell_data[cell]['peaks'][pattern] if self.cell_data[cell]['active'][pattern] == 1 else 0 for pattern in range(self.nb_patterns)])

    def temporal_plot_per_cell(self, cell):
        plt.figure()
        for pattern in range(len(self.cell_data[cell]['intensities'])):
            if self.cell_data[cell]['active'][pattern] == 1:
                plt.plot(self.cell_data[cell]['intensities'][pattern], label='Pattern ' + str(pattern) + ' (active)', linewidth=3)
                # plt.plot(self.cell_data[cell]['peak_locations'][pattern], self.cell_data[cell]['peaks'][pattern], marker='*', label='Peak of ' + str(pattern))
            else:
                plt.plot(self.cell_data[cell]['intensities'][pattern], label='Pattern ' + str(pattern) + ' (inactive)', alpha=0.2, linewidth=3)
        plt.hlines(self.cell_data[cell]['mean'] + 3*self.cell_data[cell]['std'], 0, self.nb_pattern_timepoints-1, label='Activity threshold (µ + 3σ)')
        plt.xticks(range(0,self.nb_pattern_timepoints+1,9), [str(ticklabel) for ticklabel in range(0,9)],fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (s)',fontsize=20)
        plt.ylabel('Fluorescence (dF/F)',fontsize=20)
        plt.legend(fontsize=14)
        plt.title('Temporal activity of cell ' + str(cell),fontsize=20)
        plt.tight_layout()
        # plt.show()

    def temporal_plot_per_pattern(self, pattern=0):
        if pattern > self.nb_patterns-1:
            print('Choose a pattern between 0 and ' + str(self.nb_patterns-1))
            return
        plt.figure()
        for cell in self.cell_data.keys():
            if self.cell_data[cell]['active'][pattern] == 1:
                plt.plot(self.cell_data[cell]['intensities'][pattern], label=str(pattern) + ' (active)', linewidth=3)
                plt.plot(self.cell_data[cell]['peak_locations'][pattern], self.cell_data[cell]['peaks'][pattern], marker='*')
            else:
                plt.plot(self.cell_data[cell]['intensities'][pattern], label=str(pattern) + '(inactive)', alpha=0.1, linewidth=3)
        plt.hlines(self.cell_data[cell]['mean'] + 3*self.cell_data[cell]['std'], 0, self.nb_pattern_timepoints-1, label='Activity threshold (µ + 3σ)')
        plt.xticks(range(0,self.nb_pattern_timepoints+1,9), [str(ticklabel) for ticklabel in range(0,9)],fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (s)',fontsize=20)
        plt.ylabel('Fluorescence (dF/F)',fontsize=20)
        plt.legend(fontsize=14)
        plt.title('Temporal activity of all cells for pattern ' + str(pattern),fontsize=20)
        plt.tight_layout()
        # plt.show()

    def temporal_plot_full_traces(self, cells = [1, 2, 3, 4, 5]):
        fig, axs = plt.subplots(len(cells), sharex=True)
        fig.subplots_adjust(hspace=0)
        ax_counter = 0
        for cell in cells:
            axs[ax_counter].plot([sample for trace in self.cell_data[cell]['intensities'] for sample in trace])
            axs[ax_counter].hlines(self.cell_data[cell]['mean'] + 3*self.cell_data[cell]['std'], 0, self.nb_pattern_timepoints*self.nb_patterns, linestyles='dashed')
            for pattern in range(self.nb_patterns):
                axs[ax_counter].vlines(pattern*self.nb_pattern_timepoints, 0, 1, color='red')
            ax_counter += 1

        fig.suptitle("Fluorescence (dF/F) of different cells")
        # plt.show()

    def correlation(self, pattern1=0, pattern2=1, plot=True):
        '''
        Use the Wilcoxon signed-rank test to get a p-value as index of separability between the two neuronal populations.
        '''

        if pattern1 == pattern2:
            print('Choose 2 different patterns')
            return

        # todo: now fluorescence of non-active cell is used, very similar results should be seen when 0 is used

        # fluorescence of non-active cell is used now:
        #fluorescence1 = [self.cell_data[cell]['peaks'][pattern1] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]
        #fluorescence2 = [self.cell_data[cell]['peaks'][pattern2] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]
        # value 0 is used for non-active cell:
        fluorescence1 = [self.cell_data[cell]['peaks'][pattern1] if self.cell_data[cell]['active'][pattern1] == 1 else 0 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or (self.cell_data[cell]['active'][pattern1] == 0 and self.cell_data[cell]['active'][pattern2] == 1)]
        fluorescence2 = [self.cell_data[cell]['peaks'][pattern2] if self.cell_data[cell]['active'][pattern2] == 1 else 0 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern2] == 1 or (self.cell_data[cell]['active'][pattern2] == 0 and self.cell_data[cell]['active'][pattern1] == 1)]

        color_value_X = [self.cell_data[cell]['centroid'][0]/397 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]
        color_value_Y = [self.cell_data[cell]['centroid'][1]/380 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]


        # signed_rank = wilcoxon(fluorescence1, fluorescence2) # Apply Wilcoxon test
        # corr_coeff = pearsonr(fluorescence1, fluorescence2, alternative='two-sided')
        corr_coeff = spearmanr(fluorescence1, fluorescence2, alternative='two-sided')
        # print('P-value for Wilcoxon signed-rank test for stim patterns ' + str(pattern1) + ' and ' + str(pattern2) + ' is ' + str(round(signed_rank.pvalue,8)))


        if plot == True:
            color_value='blue' # If no gradient should be plotted. Comment out if gradient should be plotted

            plt.figure()
            plt.axline(xy1=(0,0), slope=1, linewidth=2, color='blue')
            plt.scatter(fluorescence1, fluorescence2, c=color_value, s=30, cmap='viridis')
            xlim = plt.xlim()
            ylim = plt.ylim()
            mintick = min([xlim[0], ylim[0]])
            maxtick = max([xlim[1], ylim[1]])
            plt.xlim([mintick, maxtick])
            plt.ylim([mintick, maxtick])
            plt.text(mintick + 7/8*mintick/maxtick, mintick+0.02, 'r = ' + str(round(corr_coeff.statistic,2)), fontsize=15)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('Pattern ' + str(pattern1) + ' (dF/F)', fontsize=20)
            plt.ylabel('Pattern ' + str(pattern2) + ' (dF/F)', fontsize=20)
            plt.gca().set_aspect('equal')
            plt.title('Correlation of fluorescence', fontsize=20)
            plt.tight_layout()

        return corr_coeff.statistic

    def overlap(self, pattern1=0, pattern2=1, plot=True):
        '''
        Use the Wilcoxon signed-rank test to get a p-value as index of separability between the two neuronal populations.
        '''

        if pattern1 == pattern2:
            print('Choose 2 different patterns')
            return

        # 1 if active, 0 if inactive but active for other pattern, ignore if inactive for both
        fluorescence1 = [1 if self.cell_data[cell]['active'][pattern1] == 1 else 0 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or (self.cell_data[cell]['active'][pattern1] == 0 and self.cell_data[cell]['active'][pattern2] == 1)]
        fluorescence2 = [1 if self.cell_data[cell]['active'][pattern2] == 1 else 0 for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern2] == 1 or (self.cell_data[cell]['active'][pattern2] == 0 and self.cell_data[cell]['active'][pattern1] == 1)]
        overlap = [1 if (fluorescence1[i] == 1 and fluorescence2[i] == 1) else 0 for i in range(len(fluorescence1))]

        active_neurons_total = len(fluorescence1)
        active_neurons_1 = fluorescence1.count(1)
        active_neurons_2 = fluorescence2.count(1)
        active_neurons_overlap = overlap.count(1)

        if active_neurons_total > 0:
            return np.round(active_neurons_overlap/active_neurons_total,2), active_neurons_total
        else:
            return 0, 0

    def bias(self, pattern1=0, pattern2=1, plot=True):
        '''
        Use the Wilcoxon signed-rank test to get a p-value as index of separability between the two neuronal populations.
        '''

        if pattern1 == pattern2:
            print('Choose 2 different patterns')
            return

        differences = [self.cell_data[cell]['peaks'][pattern2] - self.cell_data[cell]['peaks'][pattern1] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]
        # normalized_differences = [(difference - min(differences))/(max(differences) - min(differences)) for difference in differences]

        centroids = [self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern1] == 1 or self.cell_data[cell]['active'][pattern2] == 1]
        x = [centroid[0] for centroid in centroids]
        y = [centroid[1] for centroid in centroids]

        if plot == True:
            plt.figure()
            plt.scatter(x, y, s=90, c=differences, cmap='viridis')
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('X Coordinate', fontsize=20)
            plt.ylabel('Y Coordinate', fontsize=20)
            plt.xlim([0, 397])
            plt.ylim([0, 380])
            plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
            plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
            plt.title('Neuron bias between pattern ' + str(pattern1) + ' and pattern ' + str(pattern2), fontsize=20)
            plt.tight_layout()

        return differences, centroids

    def spatial_analysis_bias(self, pattern1=0, pattern2=1, plot=True):
        # todo not finished yet, the subtraction from self.bias() gives a too flat spectrum of projected points, so no clear Gaussian visible anymore
        return_electrode = self.return_electrode_layer

        # Pattern 1
        pattern = pattern1
        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1]) # todo: take inactive cells into account? If yes, with 0 or with their baseline fluorescence?
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1])

        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=self.central_electrode, point2=return_electrode, plot=plot)
        if self.get_nb_active_neurons(pattern=pattern1)[0] < 20 or self.get_nb_active_neurons(pattern=pattern2)[0] < 20:
            print('Attention: only ' + str(self.get_nb_active_neurons(pattern=pattern)) + ' active neurons in pattern ' + str(pattern) + ', too low for Gaussian fit. Changing to stdev.')
            centroid_along_layer1, stdev_along_layer1 = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot)
        else:
            [centroid_along_layer1, stdev_along_layer1, gaussian_amp, gaussian_offset] = self.fit_neurons_gaussian(projected_coordinates, fluorescence, plot=plot)

        # Pattern 2
        pattern = pattern2
        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1]) # todo: take inactive cells into account? If yes, with 0 or with their baseline fluorescence?
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1])

        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=self.central_electrode, point2=return_electrode, plot=plot)
        if self.get_nb_active_neurons(pattern=pattern1)[0] < 20 or self.get_nb_active_neurons(pattern=pattern2)[0] < 20:
            print('Attention: only ' + str(self.get_nb_active_neurons(pattern=pattern)) + ' active neurons in pattern ' + str(pattern) + ', too low for Gaussian fit. Changing to stdev.')
            centroid_along_layer2, stdev_along_layer2 = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot)
        else:
            [centroid_along_layer2, stdev_along_layer2, gaussian_amp, gaussian_offset] = self.fit_neurons_gaussian(projected_coordinates, fluorescence, plot=plot)

        # Differences
        fluorescence, coordinates = self.bias(pattern1=pattern1, pattern2=pattern2, plot=plot)
        fluorescence = [value/2+0.5 for value in fluorescence] # Normalize fluorescence values to 0 and 1
        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=self.central_electrode, point2=return_electrode, plot=plot)
        if self.get_nb_active_neurons(pattern=pattern1)[0] < 20 or self.get_nb_active_neurons(pattern=pattern2)[0] < 20:
            print('Attention: only ' + str(self.get_nb_active_neurons(pattern=pattern)) + ' active neurons in pattern ' + str(pattern) + ', too low for Gaussian fit. Changing to stdev.')
            centroid_along_layer3, stdev_along_layer3 = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot)
        else:
            [centroid_along_layer3, stdev_along_layer3, gaussian_amp, gaussian_offset] = self.fit_neurons_gaussian(projected_coordinates, fluorescence, plot=plot)

        print(centroid_along_layer1, centroid_along_layer2, centroid_along_layer3)
        return centroid_along_layer1, centroid_along_layer2, centroid_along_layer3


    def correlation_per_angle(self, patterns=[0,1,2]):
        angles = []
        correlations = []
        overlaps = []
        for i in range(len(patterns)):
            for j in range(i+1,len(patterns)):
                angles.append(self.get_electrode_angles(self.central_electrode, self.return_electrodes[patterns[i]], self.return_electrodes[patterns[j]]))
                correlations.append(self.correlation(pattern1=patterns[i], pattern2=patterns[j], plot=False))
                overlaps.append(self.overlap(pattern1=patterns[i], pattern2=patterns[j]))

        angles, correlations, overlaps = zip(*sorted(zip(angles, correlations, overlaps))) # Sort the correlations in ascending order
        return angles, correlations, overlaps

    @staticmethod
    def get_electrode_angles(central, elec1, elec2):
        P12 = np.sqrt((central[0]-elec1[0])**2 + (central[1]-elec1[1])**2)
        P13 = np.sqrt((central[0]-elec2[0])**2 + (central[1]-elec2[1])**2)
        P23 = np.sqrt((elec1[0]-elec2[0])**2 + (elec1[1]-elec2[1])**2)
        return np.rad2deg(np.arccos((P12**2 + P13**2 - P23**2)/(2*P12*P13)))

    def plot_activity(self, pattern=0, neurons='all', arrow=False):
        plt.figure()
        for cell in self.cell_data.keys():
            if neurons == 'all':
                cell_intensity = self.cell_data[cell]['peaks'][pattern] if self.cell_data[cell]['active'][pattern] == True else 0
            elif neurons == 'unique':
                cell_intensity = self.cell_data[cell]['peaks'][pattern] - self.cell_data[cell]['common_activity'] if self.cell_data[cell]['active'][pattern] == True else 0
            elif neurons == 'common':
                cell_intensity = self.cell_data[cell]['common_activity']
            plt.scatter(self.cell_data[cell]['centroid'][0], self.cell_data[cell]['centroid'][1], s=90, c="blue", alpha=cell_intensity/self.max_activity)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('X Coordinate', fontsize=20)
        plt.ylabel('Y Coordinate', fontsize=20)
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        if arrow==True:
            plt.arrow(self.central_electrode[0], self.central_electrode[1], self.return_electrodes[pattern][0]-self.central_electrode[0], self.return_electrodes[pattern][1]-self.central_electrode[1], length_includes_head=True, width=5, head_width=20, color='black')
        plt.title('Active neurons (' + neurons + ') for pattern ' + str(pattern), fontsize=20)
        plt.tight_layout()
        # plt.show()

    def two_dimensional_kernel_density_estimate(self, pattern=0): # Currently not in use anymore
        '''
        2D Kernel Density Estimate of the data
        '''
        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys()])
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] if self.cell_data[cell]['active'][pattern] == 1 else 0 for cell in self.cell_data.keys()])

        kde = KernelDensity(bandwidth=100, kernel='gaussian') # Choose model and parameters
        kde.fit(coordinates, sample_weight=fluorescence) # Train model

        grid_size = 100 # 100 points in x and in y direction # Set to a larger grid size for a higher-resolution image as svg
        x_grid, y_grid = np.meshgrid(np.linspace(0, 397, grid_size), np.linspace(0, 380, grid_size))
        grid_points = np.vstack([x_grid.ravel(), y_grid.ravel()]).T
        density = np.exp(kde.score_samples(grid_points)).reshape(x_grid.shape) # Evaluate model for all points on the grid

        fig = plt.figure()
        plt.pcolormesh(x_grid, y_grid, density, shading='auto')
        plt.scatter(coordinates[:,0], coordinates[:,1], c=fluorescence/self.max_activity, cmap='viridis', edgecolors='k', linewidths=1)
        plt.colorbar(label='Values')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.xlim([0, 397])
        plt.ylim([0, 380])
        plt.gca().invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        # plt.legend()
        plt.title('2D Kernel Density Estimate for stim. pattern ' + str(pattern))
        # plt.show()
        plt.close()

        return coordinates, fluorescence, x_grid, y_grid, density

    def projected_kernel_density_estimate(self, pattern=0, point1=(0,0), point2=(397,380), gaussian_fit=False): # Currently not in use anymore
        '''
        Projection of the data before the Kernel Density Estimate is useful when trying the understand the spatial
        distribution of cell activity along a specific axis. The density estimate now reflects the distribution
        of the data along the projected data.
        If projection would only take place after the Kernel Density Estimate, then you just integrate the density
        function, but the density estimate is still based on a 2D-distribution.
        '''
        direction_vector = np.array([point2[0] - point1[0], point2[1] - point1[1]])
        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1]) # todo: take inactive cells into account? If yes, with 0 or with their baseline fluorescence?
        projected_points = np.array([np.dot(np.array(centroid) - np.array(point1), direction_vector) / np.dot(direction_vector, direction_vector) for centroid in coordinates])
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1])

        # Perform kernel density estimation
        kde = KernelDensity(bandwidth=0.6, kernel='gaussian') # Grid only goes from about -2 to +2 so take smaller bandwidth than 2D-KDE where grid goes from 0 to 397
        kde.fit(projected_points.reshape(-1, 1), sample_weight=fluorescence) # Train model

        # Define grid along the projected axis
        grid_size = 100
        grid = np.linspace(min(projected_points), max(projected_points), grid_size).reshape(-1, 1)
        density = np.exp(kde.score_samples(grid))

        # Plot projected KDE

        fig, ax = plt.subplots(4, 1, figsize=(5, 15))
        ax = ax.flatten()

        for cell in range(len(coordinates)):
            projected_coordinate = [point1[0] + coordinates[cell] * direction_vector[0], point1[1] + coordinates[cell] * direction_vector[1]]
            ax[0].scatter(coordinates[cell][0], coordinates[cell][1], s=90, c="blue", alpha=fluorescence[cell]/self.max_activity)
            ax[0].scatter(projected_coordinate[0], projected_coordinate[1], s=10, c="red", alpha=fluorescence[cell]/self.max_activity)

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

        # Plot 1D density function along the projected axis on ax1
        ax[2].plot(grid, density, color='red', linestyle='-')
        ax[2].set_xlabel('Distance along the projected axis')
        ax[2].set_ylabel('Density')
        ax[2].set_title('1D Kernel Density Estimate along Projected Axis')

        ax[2].set_xticks([0,1],["Point 1", "Point 2"])
        ax[2].sharex(ax[1])

        if gaussian_fit == True:
            projected_points, fluorescence = zip(*sorted(zip(projected_points, fluorescence))) # Sort the projected_points in ascending order

            nb_bins = 30
            # Create 10 equidistant bins
            bins = np.linspace(projected_points[0], projected_points[-1], nb_bins+1)  # 11 edges for 10 bins
            cumulative_values = np.zeros(nb_bins)

            bin_index = 0
            # Accumulate the fluorescence values in the bins
            for i in range(len(projected_points)):
                # Check if the projected point exceeds the current bin edge
                while bin_index < nb_bins-1 and projected_points[i] > bins[bin_index + 1]:
                    bin_index += 1
                # Add the fluorescence value to the corresponding bin
                cumulative_values[bin_index] += fluorescence[i]
            
            # Prepare the new list with mid-points of the bins and their corresponding cumulative values
            projected_points = (bins[:-1] + bins[1:]) / 2
            fluorescence = cumulative_values

            x0_estimated = grid[np.argmax(density)][0] # Take the maximum of the KDE as initial guess for the peak of the Gaussian
            p0 = [x0_estimated, 1, max(fluorescence), min(fluorescence)] # Initial guess for p = [mu, sigma, amplitude, offset]
            p_bounds = ((x0_estimated*0.8, 0, min(fluorescence), min(fluorescence)*0.5), (x0_estimated*1.2, np.inf, max(fluorescence)*2, max(fluorescence))) # Boundaries for the different parameters
            popt, pcov = curve_fit(self.gauss, projected_points, fluorescence, p0 =p0, bounds=p_bounds, maxfev=5000, method='trf', loss='linear') # Fit the points and their values to the given function todo possibly choose different loss function for regularization against outliers (eg cauchy)
            fitted_curve = self.gauss(projected_points, popt[0], popt[1], popt[2], popt[3]) # Compute the fitted curve for the projected_points
        else:
            fitted_curve = fluorescence
            popt = [None, None, None, None]
        ax[3].scatter(projected_points, fluorescence)
        ax[3].plot(projected_points, fitted_curve)

        plt.tight_layout()
        plt.show()
        plt.close()

        print(popt[0], popt[1])

        return grid, density, popt[0], popt[1], popt[2]

    @staticmethod
    def gauss(x, x0, sigma, amp, offset):
        return offset + amp * np.exp(-(x-x0)**2/(2*sigma**2))

    def kde(self, pattern=0, plot_arrow=True, plot_axes=True, plot_ellipse=False, plot_kde=True):
        coordinates, fluorescence, x_grid, y_grid, xy_density = self.two_dimensional_kernel_density_estimate(pattern=pattern)

        # Projection along the layer: can be modelled by a Gaussian
        proj_layer_grid, proj_layer_density, mu_layer, sigma_layer, amp_layer = self.projected_kernel_density_estimate(pattern=pattern, point1=self.central_electrode, point2=self.return_electrode_layer, gaussian_fit=True)
        max_layer = proj_layer_grid[np.argmax(proj_layer_density)][0] # Compute the grid value that corresponds to the index of the highest density value
        max_layer_xy = (self.central_electrode[0] + max_layer*(self.return_electrode_layer[0]-self.central_electrode[0]), self.central_electrode[1] + max_layer*(self.return_electrode_layer[1]-self.central_electrode[1])) # Compute the X and Y for this maximal value
        mu_layer_xy = (self.central_electrode[0] + mu_layer*(self.return_electrode_layer[0]-self.central_electrode[0]), self.central_electrode[1] + mu_layer*(self.return_electrode_layer[1]-self.central_electrode[1])) # Compute the X and Y for this maximal value
        sigma_layer_xy = sigma_layer*np.sqrt((self.return_electrode_layer[0]-self.central_electrode[0])**2 + (self.return_electrode_layer[1]-self.central_electrode[1])**2) # Convert sigma from projected points to sigma in FOV

        # Projection along the column: more uniform, can't be modelled by a Gaussian, take the peak of the KDE as center of activity
        proj_col_grid, proj_col_density, _, _, _ = self.projected_kernel_density_estimate(pattern=pattern, point1=self.central_electrode, point2=self.return_electrode_column, gaussian_fit=False)
        max_col = proj_col_grid[np.argmax(proj_col_density)][0]
        max_col_xy = (self.central_electrode[0] + max_col*(self.return_electrode_column[0]-self.central_electrode[0]), self.central_electrode[1] + max_col*(self.return_electrode_column[1]-self.central_electrode[1]))

        # Find intersection of mu_layer_xy (peak of fitted Gaussian) and max_col_xy (peak of KDE)

        # Line along layer through max_col_xy: a1*y+b1*x+c1=0
        a1 = self.return_electrode_layer[1]-self.central_electrode[1]
        b1 = -(self.return_electrode_layer[0]-self.central_electrode[0])
        c1 = -b1*max_col_xy[1] - a1*max_col_xy[0]

        # Line along column through mu_layer_xy: a2*y+b2*x+c2=0
        a2 = self.return_electrode_column[1]-self.central_electrode[1]
        b2 = -(self.return_electrode_column[0]-self.central_electrode[0])
        # c2 = -b2*max_layer_xy[1] - a2*max_layer_xy[0]
        c2 = -b2*mu_layer_xy[1] - a2*mu_layer_xy[0]

        max_xy = ((b1*c2-b2*c1)/(a1*b2-a2*b1), (c1*a2-c2*a1)/(a1*b2-a2*b1))

        # plt.close('all') # Close all other figures

        if plot_kde == True:
            fig = plt.figure(figsize=(8,12))
            ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=1, rowspan=2)  # 1st row, 1st column, spanning 1 column
            ax2 = plt.subplot2grid((3, 2), (0, 1), colspan=1, rowspan=2)  # 1st row, 2nd column, spanning 1 column
            ax3 = plt.subplot2grid((3, 2), (2, 0), colspan=2)  # 2nd row, 1st column, spanning 2 columns
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1)

        if plot_axes==True:
            # ax1.plot([central_elec[0], return_elec_layer[0]],[central_elec[1], return_elec_layer[1]], 'limegreen', linestyle="--", label='Along layer IV')
            ax1.axline(self.central_electrode, self.return_electrode_layer, color='limegreen', label='Along layer')
            # ax1.plot([central_elec[0], return_elec_col[0]],[central_elec[1], return_elec_col[1]], 'darkgreen', linestyle="--", label='Along cortical column')
            ax1.axline(self.central_electrode, self.return_electrode_column, color='darkgreen', label='Along column')
        ax1.scatter(coordinates[:,0], coordinates[:,1], s=90, c="blue", alpha=fluorescence/self.max_activity) # Only plots the active neurons (in contrast to plot_activity)
        ax1.scatter(self.central_electrode[0], self.central_electrode[1], color='orange', s=110, marker='s', label='Central elec. (L4)', zorder=3)
        # ax1.scatter(return_elec_col[0], return_elec_col[1], color='darkgreen', label='Adjacent elec. (L2/3)')
        # ax1.scatter(return_elec_layer[0], return_elec_layer[1], color='limegreen', label='Adjacent elec. (L4)')
        ax1.scatter(self.return_electrodes[pattern][0], self.return_electrodes[pattern][1], color='gold', s=110, marker='s', label='Return elec.', zorder=3)
        if plot_arrow==True:
            ax1.arrow(self.central_electrode[0], self.central_electrode[1], self.return_electrodes[pattern][0]-self.central_electrode[0], self.return_electrodes[pattern][1]-self.central_electrode[1], length_includes_head=True, width=5, head_width=20, color='black', zorder=3)
            if self.central_electrode == self.return_electrodes[pattern]: # Otherwise an arrow with length 0 is plotted
                circular_arrow = FancyArrowPatch((self.central_electrode[0]-10, self.central_electrode[1]-10), (self.central_electrode[0]+10, self.central_electrode[1]+10), arrowstyle='simple,head_width=10,head_length=12', connectionstyle='arc3,rad=2.5', zorder=3, color='black', linewidth=4)
                ax1.add_artist(circular_arrow)
        # ax1.scatter(max_layer_xy[0], max_layer_xy[1], color='red', marker='*', s=120, label='Max density', zorder=3)
        # ax1.scatter(max_col_xy[0], max_col_xy[1], color='red', marker='*', s=120, zorder=3)
        ax1.scatter(max_xy[0], max_xy[1], color='red', marker='X', s=120, zorder=3, label='Center of activity')
        # ax1.scatter(mu_xy[0], mu_xy[1], color='yellow', marker='X', s=120, zorder=3, label='Center of activity')
        if plot_ellipse==True:
            ellipse = Ellipse(max_xy, sigma_layer_xy, 200, angle=np.degrees(np.arctan(b2/a2)), alpha=0.5)
            ax1.add_artist(ellipse)
        ax1.set_xlabel('X Coordinate', fontsize=20)
        ax1.set_ylabel('Y Coordinate', fontsize=20)
        ax1.tick_params(axis='x', labelsize=15)
        ax1.tick_params(axis='y', labelsize=15)
        ax1.set_xlim([0, 397])
        ax1.set_ylim([0, 380])
        ax1.invert_yaxis()  # Invert y-axis for better comparison with ImageJ
        ax1.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
        ax1.legend(fontsize='8', loc='upper right')

        if plot_kde == True:
            pcm = ax2.pcolormesh(x_grid, y_grid, xy_density, shading='auto')
            ax2.scatter(coordinates[:,0], coordinates[:,1], c=fluorescence, cmap='viridis', edgecolors='k', linewidths=1)
            fig.colorbar(pcm, ax=ax2, label='Fluorescence intensity dF/F')
            # ax2.set_xlabel('X Coordinate')
            # ax2.set_ylabel('Y Coordinate')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim([0, 397])
            ax2.set_ylim([0, 380])
            ax2.invert_yaxis()  # Invert y-axis for better comparison with ImageJ
            ax2.set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
            ax2.set_title('2D KDE', fontsize=20)

            ax3.plot(proj_layer_grid, self.gauss(proj_layer_grid, mu_layer, sigma_layer, 1, 0), color='red', linestyle='-')
            ax3.set_xlabel('Distance along layer IV')
            ax3.set_ylabel('Fluorescence')
            ax3.set_xticks([0,1],["Central elec. (L4)", "Adjacent elec. (L4)"])
            ax3.set_title('Gaussian distribution along cortical layer')

        fig.suptitle('Activity for stimulation pattern ' + str(pattern), fontsize=20)
        plt.tight_layout(h_pad=4)
        # plt.show()

        return max_xy, max_layer_xy, max_col_xy

    def spatial_analysis_polarity(self, pattern=0, first_electrode=None, second_electrode=(100,100), plot=True):
        if first_electrode==None:
            first_electrode=self.central_electrode
        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1]) # todo: take inactive cells into account? If yes, with 0 or with their baseline fluorescence?
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1])
        if len(coordinates) == 0:
            print('No active neurons so no spatial analysis possible')
            return 100 # can be filtered out afterwards

        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=first_electrode, point2=second_electrode, plot=plot) # projected coordinate 0 = cathodic-first, projected coordinate 1 = anodic-first
        centroid, stdev = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot) # Only do fit with stdev, because it's projected on different axes, Gaussian fit might not always be suitable
        return centroid

    def spatial_analysis(self, pattern=0, fitting_algorithm='gaussian', plot=True):
        # gaussian, kde, centroid

        coordinates = np.array([self.cell_data[cell]['centroid'] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1]) # todo: take inactive cells into account? If yes, with 0 or with their baseline fluorescence?
        fluorescence = np.array([self.cell_data[cell]['peaks'][pattern] for cell in self.cell_data.keys() if self.cell_data[cell]['active'][pattern] == 1])

        if len(coordinates) == 0:
            print('No active neurons so no spatial analysis possible')
            return None, None, None, None

        # Centroid along cortical column:
        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=self.central_electrode, point2=self.return_electrode_column, plot=False)
        centroid_along_column = self.fit_neurons_kde(projected_coordinates, fluorescence, plot=False)
        centroid_along_column_2d = (self.central_electrode[0] + centroid_along_column*(self.return_electrode_column[0]-self.central_electrode[0]), self.central_electrode[1] + centroid_along_column*(self.return_electrode_column[1]-self.central_electrode[1])) # Compute the X and Y for this maximal value

        # Centroid along cortical layer:

        # Project 2D neurons on a certain 1D axis
        projected_coordinates, fluorescence = self.project_neurons(coordinates, fluorescence, point1=self.central_electrode, point2=self.return_electrode_layer, plot=plot)
        # Fit projected neurons
        if fitting_algorithm == 'gaussian':
            if self.get_nb_active_neurons(pattern=pattern)[0] < 20:
                print('Attention: only ' + str(self.get_nb_active_neurons(pattern=pattern)) + ' active neurons in pattern ' + str(pattern) + ', too low for Gaussian fit. Changing to stdev.')
                centroid_along_layer, stdev_along_layer = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot)
            else:
                [centroid_along_layer, stdev_along_layer, gaussian_amp, gaussian_offset] = self.fit_neurons_gaussian(projected_coordinates, fluorescence, plot=plot)
        elif fitting_algorithm == 'kde':
            centroid_along_layer = self.fit_neurons_kde(projected_coordinates, fluorescence, plot=plot)
            stdev_along_layer = None
        elif fitting_algorithm == 'stdev':
            centroid_along_layer, stdev_along_layer = self.fit_neurons_stdev(projected_coordinates, fluorescence, plot=plot)
        else:
            print('Choose an existing fitting algorithm: gaussian/kde/stdev')
            return

        centroid_along_layer_2d = (self.central_electrode[0] + centroid_along_layer*(self.return_electrode_layer[0]-self.central_electrode[0]), self.central_electrode[1] + centroid_along_layer*(self.return_electrode_layer[1]-self.central_electrode[1])) # Compute the X and Y for this maximal value
        # Project center on projection axis back to 2D to get 2D-centroid
        centroid_2d = self.intersect_projections(centroid_along_column_2d, centroid_along_layer_2d)

        return centroid_along_column, centroid_along_layer, stdev_along_layer, centroid_2d

    def project_neurons(self, original_points, amplitude, point1=(0,0), point2=(397,380), plot=False):
        direction_vector = np.array([point2[0] - point1[0], point2[1] - point1[1]])
        projected_points = np.array([np.dot(np.array(centroid) - np.array(point1), direction_vector) / np.dot(direction_vector, direction_vector) for centroid in original_points])

        if len(original_points) == 0:
            print('No active neurons so no projection possible.')
            return

        projected_points, amplitude = zip(*sorted(zip(projected_points, amplitude))) # Sort the projected_points in ascending order
        projected_points = np.array(projected_points)
        amplitude = np.array(amplitude)

        if plot==True:
            fig, ax = plt.subplots(2, 1)
            ax = ax.flatten()

            # Plot for original points, line, and projected points
            ax[0].plot([point1[0], point2[0]], [point1[1], point2[1]], color='blue', label='Line')
            for cell in range(len(original_points)):
                projected_coordinate = [point1[0] + projected_points[cell] * direction_vector[0], point1[1] + projected_points[cell] * direction_vector[1]]
                ax[0].scatter(original_points[cell][0], original_points[cell][1], s=90, c="blue", alpha=amplitude[cell]/self.max_activity)
                ax[0].scatter(projected_coordinate[0], projected_coordinate[1], s=20, c="red", alpha=amplitude[cell]/self.max_activity)

            ax[0].set_xlabel('X Coordinate')
            ax[0].set_ylabel('Y Coordinate')
            ax[0].set_xlim([0, 397])
            ax[0].set_ylim([0, 380])
            ax[0].invert_yaxis()  # Invert y-axis for better comparison with ImageJ
            ax[0].set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
            ax[0].set_title('Projection of points onto projection axis')
            ax[0].legend()

            # Plot for points along the projection axis
            ax[1].scatter(projected_points, amplitude, color='blue')
            ax[1].set_xticks([])
            ax[1].set_ylabel('Amplitude')
            ax[1].set_title('Points along projection axis')
            ax[1].set_aspect('auto')

            fig.tight_layout()

        return projected_points, amplitude

    def fit_neurons_gaussian(self, points, amplitude, plot=False):
        # todo only with neurons?
        x0_estimated, mu0_estimated = self.fit_neurons_stdev(points, amplitude, plot=False) # Take the mean and standard deviation as initial guess for the peak and standard deviation of the Gaussian

        nb_bins = 20 # or set to a fixed size
        # Create 10 equidistant bins
        bins = np.linspace(points[0], points[-1], nb_bins+1)  # 11 edges for 10 bins
        cumulative_values = np.zeros(nb_bins)

        bin_index = 0
        # Accumulate the fluorescence values in the bins
        for i in range(len(points)):
            # Check if the projected point exceeds the current bin edge
            while bin_index < nb_bins-1 and points[i] > bins[bin_index + 1]:
                bin_index += 1
            # Add the fluorescence value to the corresponding bin
            cumulative_values[bin_index] += amplitude[i]

        # Prepare the new list with mid-points of the bins and their corresponding cumulative values
        points = (bins[:-1] + bins[1:]) / 2
        amplitude = cumulative_values


        p0 = [x0_estimated, mu0_estimated, max(amplitude), min(amplitude)] # Initial guess for p = [mu, sigma, amplitude, offset]
        p_bounds = ((x0_estimated*0.8 if x0_estimated >= 0 else x0_estimated*1.2, 0, min(amplitude), min(amplitude)*0.5), (x0_estimated*1.2 if x0_estimated >= 0 else x0_estimated*0.8, np.inf, max(amplitude)*2, max(amplitude))) # Boundaries for the different parameters
        popt, pcov = curve_fit(self.gauss, points, amplitude, p0 =p0, bounds=p_bounds, maxfev=5000, method='trf', loss='linear') # Fit the points and their values to the given function todo possibly choose different loss function for regularization against outliers (eg cauchy)
        fitted_curve = self.gauss(points, popt[0], popt[1], popt[2], popt[3]) # Compute the fitted curve for the projected_points

        if plot==True:
            plt.figure()
            # Plot for points along the projection axis
            plt.scatter(points, amplitude, color='blue')
            plt.plot(points, fitted_curve, color='red', linestyle='-')
            plt.xlabel('Density along projected axis')
            plt.xticks([0,1],["Point 1", "Point 2"])
            plt.title('1D Gaussian fit of projected points')

        return popt

    def fit_neurons_kde(self, points, amplitude, plot=False):
        # Perform kernel density estimation
        kde = KernelDensity(bandwidth=0.6, kernel='gaussian') # Grid only goes from about -2 to +2 so take smaller bandwidth than 2D-KDE where grid goes from 0 to 397
        kde.fit(points.reshape(-1, 1), sample_weight=amplitude) # Train model

        # Define grid along the projected axis
        grid_size = 100
        grid = np.linspace(min(points), max(points), grid_size).reshape(-1, 1)
        density = np.exp(kde.score_samples(grid))
        centroid = grid[np.argmax(density)][0]

        if plot==True:
            plt.figure()
            # Plot for points along the projection axis
            plt.scatter(points, amplitude, color='blue')
            plt.plot(grid, density, color='red', linestyle='-')
            plt.scatter(centroid, 0, color='red', linewidth=10)
            plt.xlabel('Density along projected axis')
            plt.xticks([0,1],["Point 1", "Point 2"])
            plt.title('1D KDE of projected points')

        return centroid


    def fit_neurons_stdev(self, points, amplitude, plot=False):
        # todo now it takes amplitude into account, but not sure if this is the correct way of doing this?

        centroid = np.average(points, weights=amplitude)
        stdev = np.std(points) # todo no weights yet for the standard deviation, possible with: https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy

        if plot==True:
            plt.figure()
            # Plot for points along the projection axis
            plt.scatter(points, amplitude, color='blue')
            plt.scatter(centroid, 0, color='red')
            plt.plot([centroid-stdev, centroid+stdev], [0, 0], color='red')
            plt.xlabel('Density along projected axis')
            plt.xticks([0,1],["Point 1", "Point 2"])
            plt.title('Centroid and standard deviation of projected points')

        return centroid, stdev

    def intersect_projections(self, centroid_column_xy, centroid_layer_xy):
        # Line along layer through max_col_xy: a1*y+b1*x+c1=0
        a1 = self.return_electrode_layer[1]-self.central_electrode[1]
        b1 = -(self.return_electrode_layer[0]-self.central_electrode[0])
        c1 = -b1*centroid_column_xy[1] - a1*centroid_column_xy[0]

        # Line along column through mu_layer_xy: a2*y+b2*x+c2=0
        a2 = self.return_electrode_column[1]-self.central_electrode[1]
        b2 = -(self.return_electrode_column[0]-self.central_electrode[0])
        # c2 = -b2*max_layer_xy[1] - a2*max_layer_xy[0]
        c2 = -b2*centroid_layer_xy[1] - a2*centroid_layer_xy[0]

        max_xy = ((b1*c2-b2*c1)/(a1*b2-a2*b1), (c1*a2-c2*a1)/(a1*b2-a2*b1))

        return max_xy

    def get_nb_active_neurons(self, pattern=0):
        nb_active_neurons = 0
        cumulative_activity = 0
        for cell in self.cell_data.keys():
            if self.cell_data[cell]['active'][pattern] == 1:
                nb_active_neurons += 1
                cumulative_activity += self.cell_data[cell]['peaks'][pattern]
        return nb_active_neurons, cumulative_activity



    ###################################################

    def amplitude_cumulative_activity(self, amp_range = 'all', plot=True):
        all_amplitudes = [6.6, 12.5, 15.4, 22.7, 25.6, 28.6, 31.5, 34.4]
        if amp_range == 'all':
            amplitudes = all_amplitudes
        elif amp_range == 'low':
            amplitudes = all_amplitudes[:4]
        elif amp_range == 'high':
            amplitudes = all_amplitudes[-4:]

        step = 2 # 2 for old programming (filter out asymmetric ones), 1 for new programming (only symmetric ones)
        nb_active_neurons = [0]*int(self.nb_patterns/step)
        cumulative_activity = [0]*int(self.nb_patterns/step)
        for cell in self.cell_data.keys():
            for pattern in range(len(nb_active_neurons)):
                if self.cell_data[cell]['active'][pattern*step] == 1:
                    nb_active_neurons[pattern] += 1
                    cumulative_activity[pattern] += self.cell_data[cell]['peaks'][pattern*step]

        # todo: in absolute values or in percentage?

        if plot==True:
            fig, ax1 = plt.subplots()
            ax1.plot(amplitudes, nb_active_neurons, color='blue', linewidth=3)
            ax1.set_ylabel('Nb. of active neurons', fontsize=20)
            ax1.set_xlabel('Current amplitude [µA]', fontsize=20)
            ax1.tick_params(axis='both', labelsize=15)
            ax1.set_xticks([0, 5, 10, 15, 20, 25, 30, 35, 40])
            ax1.tick_params(axis='x', labelsize=15)
            #ax2 = ax1.twinx()
            #ax2.plot(amplitudes, cumulative_activity, color='red', linewidth=3)
            #ax2.set_ylabel('Cumulative activity', color='red', fontsize=20)
            #ax2.tick_params(axis='y', labelsize=15, labelcolor='red')
            plt.title('Activity for increasing stimulation amplitude', fontsize=20)
            plt.tight_layout()

        return amplitudes, nb_active_neurons, cumulative_activity

def get_all_electrode_locations(central, return_col, return_layer):
    # Assumes a rectangular grid, always take the middle of the needle
    all_electrodes = [return_col, return_layer]
    difference_row = tuple(np.subtract(return_layer, central))
    all_electrodes.append(add_tuple(return_col, difference_row, 1/2))
    all_electrodes.append(add_tuple(central, difference_row, 2))
    all_electrodes.append(add_tuple(return_col, difference_row, 2))
    all_electrodes.append(add_tuple(return_col, difference_row, 3/2))
    all_electrodes.append(add_tuple(return_col, difference_row, 1))
    all_electrodes.append(add_tuple(return_col, difference_row, -1/2))
    all_electrodes.append(add_tuple(return_col, difference_row, -1))
    all_electrodes.append(add_tuple(return_col, difference_row, -3/2))
    all_electrodes.append(add_tuple(return_col, difference_row, -2))
    all_electrodes.append(add_tuple(central, difference_row, -2))
    all_electrodes.append(add_tuple(central, difference_row, -1))
    print('All electrodes: ', all_electrodes)
    return all_electrodes

def add_tuple(tuple1, tuple2, multiplication):
    return (int(tuple1[0] + multiplication*tuple2[0]), int(tuple1[1] + multiplication*tuple2[1]))

# Load electrodes: [central_electrode, return_electrodes, return_electrode_column, return_electrode_layer]

'''
# To load in data and check whether it's correct
date = "04-24"
name = "Slice7 - local"
patterns = 6

###
directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-" + date
experiment = ExtractFluorescence(directory=directory, recording_name=name, nb_patterns=patterns)
electrode_locations= [(241, 155), [(328, 254)]*patterns, (142, 224), (328, 254)]
comparison = ProcessFluorescence(experiment.cell_data, electrode_locations)


for pattern in range(comparison.nb_patterns):
    comparison.plot_activity(neurons='all', arrow=False, pattern=pattern)

for cell in [5, 10, 15, 20, 25, 30]:
    comparison.temporal_plot_per_cell(cell)
'''

############################################################
##################### DIRECTIONS ###########################
############################################################

list_directions20 = []
list_directions10 = []

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-26"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 1+2+3", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
jan26_2_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(jan26_2_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-26"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 4+5+6", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
jan26_2_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(jan26_2_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
jan31_4_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(jan31_4_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 3+4+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
jan31_4_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(jan31_4_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-05"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 0+4+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb05_4_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb05_4_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-05"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb05_7_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(feb05_7_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-05"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 3+4+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb05_7_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb05_7_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-06"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb06_2_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb06_2_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-06"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice5 - 0+2+4", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb06_5_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(feb06_5_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-06"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice5 - 1+3+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb06_5_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb06_5_20)

# Note: 30um above regular plane
directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-07"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 3+5+7", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb07_4_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb07_4_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice1 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb14_1_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb14_1_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb14_2_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb14_2_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 3+4+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb14_2_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(feb14_2_10)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice3 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb14_3_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb14_3_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice3 - 3+4+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb14_3_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(feb14_3_10)

# todo Possibly add feb 15: Slice2, Slice3

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-15"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb15_7_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb15_7_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-15"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice8 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb15_8_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb15_8_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-16"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice3 - 0+1+2", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
feb16_3_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(feb16_3_20)

# todo Add 23 apr slice 3 but only few directions

# todo Add 23 apr slice 7 but only few directions

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 0+7+5", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
apr24_7_20 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions20.append(apr24_7_20)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 1+6+4", nb_patterns=13)
electrode_locations = [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
apr24_7_10 = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_directions10.append(apr24_7_10)

############################################################
##################### ASYMMETRY ############################
############################################################

list_asymm = []

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - asymm", nb_patterns=24)
electrode_locations = [(241, 155), [(145, 215), (324, 258), (200, 220)]*8, (145, 215), (324, 258)]
jan31_4_asymm = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_asymm.append(jan31_4_asymm)

# todo add 5 feb slice 4 but not many neurons (alreayd typed below)

#directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-05"
#experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - asymm", nb_patterns=24)
#electrode_locations = [(241, 155), [(145, 215), (324, 258), (200, 220)]*8, (145, 215), (324, 258)]
#feb05_4_asymm = ProcessFluorescence(experiment.cell_data, electrode_locations)
#list_asymm.append(feb05_4_asymm)


# todo add 5 feb slice 7 but not many neurons

# todo add 14 feb slice 1 but only ok for P0

# todo add 14 feb slice 3 but only ok for P1

# todo add 14 feb slice 4 but probably not ok

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-15"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - asymm", nb_patterns=24)
electrode_locations = [(241, 155), [(145, 215), (324, 258), (200, 220)]*8, (145, 215), (324, 258)]
feb15_2_asymm = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_asymm.append(feb15_2_asymm)

# todo add 16 feb slice 6 but not many neurons

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-16"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice8 - asymm", nb_patterns=24)
electrode_locations = [(241, 155), [(145, 215), (324, 258), (200, 220)]*8, (145, 215), (324, 258)]
feb16_8_asymm = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_asymm.append(feb16_8_asymm)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - asymm", nb_patterns=24)
electrode_locations = [(241, 155), [(145, 215), (324, 258), (200, 220)]*8, (145, 215), (324, 258)]
apr24_7_asymm = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_asymm.append(apr24_7_asymm)

############################################################
##################### DEPTH ################################
############################################################

list_depth = []

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - depth", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
jan31_4_depth = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_depth.append(jan31_4_depth)

# todo add 5 feb slice 4 but not many neurons visible

# todo add 14 feb slice 1 but not many neurons visible

# todo add 23 apr slice 5 but not many neurons visible

# todo add 23 apr slice 7 but not many neurons visible

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice6 - depth", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
apr24_6_depth = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_depth.append(apr24_6_depth)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - depth", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
apr24_7_depth = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_depth.append(apr24_7_depth)

# todo add 24 apr slice 8 but not many neurons visible

# todo add 2 mei slice 3 but not many neurons visible

# todo add 2 mei slice 9

############################################################
##################### LOCAL ################################
############################################################

list_local = []

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - local", nb_patterns=6)
electrode_locations = [(241, 155), [(324, 258)]*6, (145, 215), (324, 258)]
jan31_4_local = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_local.append(jan31_4_local)

# todo add 5 feb slice 4 but not many neurons visible

# todo add 14 feb slice 4

# todo add 16 feb slice 8 but there's something weird in bottom left corner

# todo add 23 apr slice 5 but not many neurons visible

# todo add 23 apr slice 7 but not many neurons visible

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice6 - local", nb_patterns=6)
electrode_locations = [(241, 155), [(324, 258)]*6, (145, 215), (324, 258)]
apr24_6_local = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_local.append(apr24_6_local)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - local", nb_patterns=6)
electrode_locations = [(241, 155), [(324, 258)]*6, (145, 215), (324, 258)]
apr24_7_local = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_local.append(apr24_7_local)

# todo add 24 apr slice 8 but not many neurons visible

# todo add 2 mei slice 3 but not many neurons visible

# todo add 2 mei slice 9

############################################################
##################### AMPLITUDE ############################
############################################################

list_amplitude = []

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-26"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice2 - 7", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
jan26_2_low = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(jan26_2_low)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 14", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
jan31_4_low = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(jan31_4_low)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-05"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 7+8", nb_patterns=16)
electrode_locations = [(241, 155), [(324, 258)]*16, (145, 215), (324, 258)]
feb05_4_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb05_4_all)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-06"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 2+3", nb_patterns=16)
electrode_locations = [(241, 155), [(324, 258)]*16, (145, 215), (324, 258)]
feb06_7_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb06_7_all)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-06"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice9 - 1", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
feb06_9_high = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb06_9_high)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-07"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice4 - 1", nb_patterns=8)
electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
feb07_4_high = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb07_4_high)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice1 - 16+17", nb_patterns=16)
electrode_locations = [(241, 155), [(324, 258)]*16, (145, 215), (324, 258)]
feb14_1_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb14_1_all)

directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-02-14"
experiment = ExtractFluorescence(directory=directory, recording_name="Slice3 - 10+11", nb_patterns=16)
electrode_locations = [(241, 155), [(324, 258)]*16, (145, 215), (324, 258)]
feb14_3_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
list_amplitude.append(feb14_3_all)

# todo add 16 feb slice 2 but not sure if useful

# todo add 16 feb slice 4 but not sure if useful

# todo add 24 apr slice 6 and 7


#directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
#experiment = ExtractFluorescence(directory=directory, recording_name="Slice6 - 4", nb_patterns=8)
#electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
#apr24_6_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
#list_amplitude.append(apr24_6_all)

#directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-04-24"
#experiment = ExtractFluorescence(directory=directory, recording_name="Slice7 - 8", nb_patterns=8)
#electrode_locations = [(241, 155), [(324, 258)]*8, (145, 215), (324, 258)]
#apr24_7_all = ProcessFluorescence(experiment.cell_data, electrode_locations)
#list_amplitude.append(apr24_7_all)


# todo add 2 mei slice 9



############################################################
##################### DATA ANALYSIS ########################
############################################################


############################################################
######### DATA ANALYSIS DIRECTIONALITY SPATIAL #############
############################################################

centroids = [[], [], [], [], [], [], [], [], [], [], [], [], []]

for slice in list_directions20:
    print(slice)
    for pattern in range(0,slice.nb_patterns):
        print('Pattern: ' + str(pattern))
        nb_active_neurons = slice.get_nb_active_neurons(pattern=pattern)[0]
        print(nb_active_neurons)
        if nb_active_neurons < 30:
            continue
        centroid_column, centroid_layer, stdev, _, = slice.spatial_analysis(pattern=pattern, fitting_algorithm='gaussian', plot=False)
        print(centroid_column)
        print(centroid_layer)
        centroids[pattern].append([centroid_column, centroid_layer, stdev])

colors = ['red', 'blue', 'orange', 'purple', 'yellow', 'purple', 'grey', 'green']
return_elec = [(0, 1), (1, 0), (0.5, 1), (2, 0), (2, 1), (1.5, 1), (1, 1), (-0.5, 1), (-1, 1)]

plt.figure()
for pattern in [0, 1, 2, 6, 7]:
    for centroid in range(len(centroids[pattern])):
        plt.plot(centroids[pattern][centroid][1], centroids[pattern][centroid][0], 'o', markersize=5, color=colors[pattern], alpha=0.5, label='Intra-slice centroids')

        # plt.plot([centroids[pattern][centroid][1]-centroids[pattern][centroid][2], centroids[pattern][centroid][1]+centroids[pattern][centroid][2]], [centroids[pattern][centroid][0], centroids[pattern][centroid][0]], '-', linewidth=3, color=colors[pattern])
    plt.plot(np.mean([centroids[pattern][centroid][1] for centroid in range(len(centroids[pattern]))]), np.mean([centroids[pattern][centroid][0] for centroid in range(len(centroids[pattern]))]), 'o', markersize=10, color=colors[pattern], label='Inter-slice centroid', zorder=3)
    plt.arrow(0, 0, return_elec[pattern][0], return_elec[pattern][1], length_includes_head=True, width=0.01, head_width=0.05, color=colors[pattern])
plt.xlabel('Distance parallel to surface [µm]', fontsize=20)
plt.ylabel('Distance perpendicular \nto surface [µm]', fontsize=20)
plt.xticks(ticks=[-1, -0.66, -0.33, 0, 0.33, 0.66, 1], labels=['-150', '-100', '-50', '0', '50', '100', '150'], fontsize=15)
plt.yticks(ticks=[-0.33, 0, 0.33, 0.66, 1, 1.33, 1.66, 2, 2.33], labels=['-50', '0', '50', '100', '150', '200', '250', '300', '350'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)
plt.gca().get_legend().legendHandles[0].set_color('black')
plt.gca().get_legend().legendHandles[1].set_color('black')

plt.text(-0.95, 0.15, 'Layer 2/3', fontsize=15)
plt.hlines(0.05, -0.95, -0.25, colors='black', linestyles='dashed')
plt.text(-0.95, -0.15, 'Layer 4', fontsize=15)

plt.title('Activity centroids for different directions', fontsize=20)
plt.tight_layout()


# todo PCA toepassen op de vectoren in de 13-dimensionale ruimte van de 13 verschillende clusters om zo een significant verschillend aantal clusters te ontdekken

# down here example of how gaussian fit can be used for directionality

'''
directory = "C:/D_Drive/OneDrive - KU Leuven/NERF-VIB/Animal experiments/Experiment 2024-01-31"
experiment10 = ExtractFluorescence(directory=directory, recording_name="Slice4 - 0+1+2", nb_patterns=13)
electrode_locations= [(241, 155), get_all_electrode_locations((241, 155),(145, 215), (324, 258)), (145, 215), (324, 258)]
comparison10 = ProcessFluorescence(experiment10.cell_data, electrode_locations)
popt0 = [0,0]
grid, density, popt0[0], popt0[1], _ = comparison10.projected_kernel_density_estimate(pattern=0, point1=(241,155), point2=(324,258), gaussian_fit=True)
popt1 = [0,0]
grid, density, popt1[0], popt1[1], _ = comparison10.projected_kernel_density_estimate(pattern=1, point1=(241,155), point2=(324,258), gaussian_fit=True)
popt2 = [0,0]
grid, density, popt2[0], popt2[1], _ = comparison10.projected_kernel_density_estimate(pattern=2, point1=(241,155), point2=(324,258), gaussian_fit=True)

plt.figure()
plt.plot(grid, comparison10.gauss(grid, popt0[0], popt0[1], 1, 0), label='Pattern 0', linewidth=3)
plt.plot(grid, comparison10.gauss(grid, popt1[0], popt1[1], 1, 0), label='Pattern 1', linewidth=3)
plt.plot(grid, comparison10.gauss(grid, popt2[0], popt2[1], 1, 0), label='Pattern 2', linewidth=3)
plt.xlabel('X Coordinate along layer', fontsize=20)
plt.ylabel('Neural activity', fontsize=20)
plt.xticks(ticks=[0, 1], labels=['Central elec.', 'Elec. next needle'], fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize='8', loc='upper right')
plt.tight_layout()
plt.show()
'''

############################################################
######### DATA ANALYSIS ASYMMETRY AND POLARITY #############
############################################################

# 3 different combinations: different polarity, different asymmetry, different polarity and asymmetry

correlations_polar_20 = []
correlations_asymm_20 = []
correlations_both_20 = []
# Formula: centroid_shift = centroid from anodic-first - centroid from cathodic-first, divided by 2
# centroid_shift for direction_pattern P2 should be multiplied with ( = 150/168). Distance to virtual electrode of P2 is 168um instead of the 150um distance between the other electrodes. So to be able to compare the results of P2 with P0 and P1, a factor is required.


for slice in list_asymm:
    for direction_pattern in [1, 9, 17]: # cathodic-first and symmetric, 20uA
        if direction_pattern == 1: # P0
            second_electrode = slice.return_electrode_column
            distance_factor = 1
        elif direction_pattern == 9: # P1
            second_electrode = slice.return_electrode_layer
            distance_factor = 1
        elif direction_pattern == 17: # P2
            second_electrode = slice.return_electrodes[2]
            distance_factor = 1.12

        correlations_polar_20.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+2, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+2, plot=False),
                                      distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+2, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+0, second_electrode=second_electrode, plot=False))/2])
        correlations_polar_20.append([slice.correlation(pattern1=direction_pattern+4, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+4, pattern2=direction_pattern+6, plot=False),
                                      distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+6, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+4, second_electrode=second_electrode, plot=False))/2])

        correlations_asymm_20.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+4, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+4, plot=False)])
        correlations_asymm_20.append([slice.correlation(pattern1=direction_pattern+2, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+2, pattern2=direction_pattern+6, plot=False)])

        correlations_both_20.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+6, plot=False),
                                     distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+6, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+0, second_electrode=second_electrode, plot=False))/2])
        correlations_both_20.append([slice.correlation(pattern1=direction_pattern+2, pattern2=direction_pattern+4, plot=False), slice.overlap(pattern1=direction_pattern+2, pattern2=direction_pattern+4, plot=False),
                                     distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+2, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+4, second_electrode=second_electrode, plot=False))/2])

correlations_polar_10 = []
correlations_asymm_10 = []
correlations_both_10 = []

for slice in list_asymm:
    for direction_pattern in [0, 8, 16]: # cathodic-first and symmetric, 10uA
        if direction_pattern == 1: # P0
            second_electrode = slice.return_electrode_column
            distance_factor = 1
        elif direction_pattern == 9: # P1
            second_electrode = slice.return_electrode_layer
            distance_factor = 1
        elif direction_pattern == 17: # P2
            second_electrode = slice.return_electrodes[2]
            distance_factor = 1.12

        correlations_polar_10.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+2, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+2, plot=False),
                                      distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+2, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+0, second_electrode=second_electrode, plot=False))/2])
        correlations_polar_10.append([slice.correlation(pattern1=direction_pattern+4, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+4, pattern2=direction_pattern+6, plot=False),
                                      distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+6, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+4, second_electrode=second_electrode, plot=False))/2])

        correlations_asymm_10.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+4, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+4, plot=False)])
        correlations_asymm_10.append([slice.correlation(pattern1=direction_pattern+2, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+2, pattern2=direction_pattern+6, plot=False)])

        correlations_both_10.append([slice.correlation(pattern1=direction_pattern+0, pattern2=direction_pattern+6, plot=False), slice.overlap(pattern1=direction_pattern+0, pattern2=direction_pattern+6, plot=False),
                                     distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+6, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+0, second_electrode=second_electrode, plot=False))/2])
        correlations_both_10.append([slice.correlation(pattern1=direction_pattern+2, pattern2=direction_pattern+4, plot=False), slice.overlap(pattern1=direction_pattern+2, pattern2=direction_pattern+4, plot=False),
                                     distance_factor*(slice.spatial_analysis_polarity(pattern=direction_pattern+2, second_electrode=second_electrode, plot=False) - slice.spatial_analysis_polarity(pattern=direction_pattern+4, second_electrode=second_electrode, plot=False))/2])



threshold = 50


fig, ax = plt.subplots()
markersize=5
colors = ['blue', 'orange']
color = 'grey'
alpha = 0.5

correlations = [correlation[0] for correlation in correlations_polar_20 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(0.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha, label='20 µA')
plt.plot(0.9, np.mean(correlations), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(0.9, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[0])

correlations = [correlation[0] for correlation in correlations_polar_10 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(1.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha, label='10 µA')
plt.plot(1.1, np.mean(correlations), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(1.1, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[1])




correlations = [correlation[0] for correlation in correlations_asymm_20 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(1.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha)
plt.plot(1.9, np.mean(correlations), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(1.9, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[0])

correlations = [correlation[0] for correlation in correlations_asymm_10 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(2.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha)
plt.plot(2.1, np.mean(correlations), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(2.1, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[1])




correlations = [correlation[0] for correlation in correlations_both_20 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(2.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha)
plt.plot(2.9, np.mean(correlations), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(2.9, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[0])

correlations = [correlation[0] for correlation in correlations_both_10 if correlation[1][1] > threshold]
for dot in correlations:
    plt.plot(3.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha)
plt.plot(3.1, np.mean(correlations), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(3.1, np.mean(correlations), yerr=np.std(correlations), capsize=20, capthick=2, color=colors[1])

plt.ylim([0, 1])
plt.yticks(fontsize=15)
plt.ylabel('Spearman correlation coeff.', fontsize=20)
plt.xticks([1,2,3], ['Polarity', 'Symmetry', 'Both'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Correlation for different parameters', fontsize=20)
fig.tight_layout()



fig, ax = plt.subplots()
markersize=5
colors = ['blue', 'orange']
color = 'grey'
alpha = 0.5

overlap = [correlation[1][0] for correlation in correlations_polar_20 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(0.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha, label='20 µA')
plt.plot(0.9, np.mean(overlap), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(0.9, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[0])

overlap = [correlation[1][0] for correlation in correlations_polar_10 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(1.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha, label='10 µA')
plt.plot(1.1, np.mean(overlap), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(1.1, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[1])



overlap = [correlation[1][0] for correlation in correlations_asymm_20 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(1.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha)
plt.plot(1.9, np.mean(overlap), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(1.9, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[0])

overlap = [correlation[1][0] for correlation in correlations_asymm_10 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(2.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha)
plt.plot(2.1, np.mean(overlap), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(2.1, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[1])



overlap = [correlation[1][0] for correlation in correlations_both_20 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(2.9, dot, 'o', markersize=markersize, color=colors[0], alpha=alpha)
plt.plot(2.9, np.mean(overlap), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(2.9, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[0])

overlap = [correlation[1][0] for correlation in correlations_both_10 if correlation[1][1] > threshold]
for dot in overlap:
    plt.plot(3.1, dot, 'o', markersize=markersize, color=colors[1], alpha=alpha)
plt.plot(3.1, np.mean(overlap), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(3.1, np.mean(overlap), yerr=np.std(overlap), capsize=20, capthick=2, color=colors[1])


plt.ylim([0, 1])
plt.yticks(fontsize=15)
plt.ylabel('Binary overlap', fontsize=20)
plt.xticks([1,2,3], ['Polarity', 'Symmetry', 'Both'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Neural overlap for different parameters', fontsize=20)
fig.tight_layout()



fig, ax = plt.subplots()
markersize=5
colors = ['blue', 'orange']
alpha = 0.5

bias = [correlation[2] for correlation in correlations_polar_20 if correlation[1][1] > threshold] + [correlation[2] for correlation in correlations_both_20 if correlation[1][1] > threshold]
for centroid in bias:
    plt.plot(1, centroid, 'o', markersize=markersize, color=colors[0], alpha=alpha, label='20 µA')
plt.plot(1, np.mean(bias), 'o', markersize=markersize*2, color=colors[0])
ax.errorbar(1, np.mean(bias), yerr=np.std(bias), capsize=20, capthick=2, color=colors[0])

bias = [correlation[2] for correlation in correlations_polar_10 if correlation[1][1] > threshold] + [correlation[2] for correlation in correlations_both_10 if correlation[1][1] > threshold]
for centroid in bias:
    plt.plot(2, centroid, 'o', markersize=markersize, color=colors[1], alpha=alpha, label='10 µA')
plt.plot(2, np.mean(bias), 'o', markersize=markersize*2, color=colors[1])
ax.errorbar(2, np.mean(bias), yerr=np.std(bias), capsize=20, capthick=2, color=colors[1])

#plt.ylim([-0.1, 1.1])
plt.yticks([-0.13, 0, 0.13, 0.27, 0.4], ['-20', '0', '20', '40', '60'], fontsize=15)
plt.ylabel('Bias to cath.-first elec. [µm]', fontsize=20)
plt.xticks([1,2], ['20 µA', '10 µA'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Centroid shift towards cathodic-first ', fontsize=20)
plt.tight_layout()


############################################################
######### DATA ANALYSIS DIRECTIONALITY CELLULAR ############
############################################################

# overlap and correlation only possible (in first instance) for patterns within the same recording

threshold_10 = 80
threshold_20 = 100

'''
patterns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
plt.figure()
for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold:
            plt.plot(angles[combination], correlations[combination], '.', markersize=10, color='blue')

plt.figure()
for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold:
            plt.plot(angles[combination], overlaps[combination][0], '.', markersize=10, color='blue')
'''

markersize = 7
colors = ['blue', 'orange']

plt.figure()

patterns = [0, 1, 3, 9, 12]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[0], label='20 µA')

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[1], label='10 µA')

patterns = [5, 6, 7, 10]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[0])

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[1])

patterns = [2, 4, 8, 11]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[0])

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], correlations[combination], 'o', markersize=markersize, color=colors[1])

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.xlabel('Angle between 2 directions [°]', fontsize=20)
plt.ylabel('Spearman correlation coeff. \nof 2 directions', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
bottom, top = plt.ylim()
plt.ylim(bottom, 1)
plt.title('Correlation for different directions', fontsize=20)
plt.tight_layout()

plt.figure()

patterns = [0, 1, 3, 9, 12]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[0], label='20 µA')

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[1], label='10 µA')

patterns = [5, 6, 7, 10]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[0])

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[1])

patterns = [2, 4, 8, 11]

for slice in list_directions20:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_20:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[0])

for slice in list_directions10:
    angles, correlations, overlaps = slice.correlation_per_angle(patterns=patterns)
    for combination in range(len(angles)):
        if overlaps[combination][1] > threshold_10:
            plt.plot(angles[combination], overlaps[combination][0], 'o', markersize=markersize, color=colors[1])

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.xlabel('Angle between 2 directions [°]', fontsize=20)
plt.ylabel('Binary overlap of 2 directions', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
bottom, top = plt.ylim()
plt.ylim(bottom, 1)
plt.title('Neural overlap for different directions', fontsize=20)
plt.tight_layout()


############################################################
################ DATA ANALYSIS AMPLITUDE ###################
############################################################

plt.figure()
linewidth=2.5
amp, neurons, _ = list_amplitude[0].amplitude_cumulative_activity(amp_range='low', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
#amp, neurons, _ = list_amplitude[1].amplitude_cumulative_activity(amp_range='low', plot=False)
#plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[2].amplitude_cumulative_activity(amp_range='all', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[3].amplitude_cumulative_activity(amp_range='all', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[4].amplitude_cumulative_activity(amp_range='high', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[5].amplitude_cumulative_activity(amp_range='high', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[6].amplitude_cumulative_activity(amp_range='all', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')
amp, neurons, _ = list_amplitude[7].amplitude_cumulative_activity(amp_range='all', plot=False)
plt.plot(amp, neurons, linewidth=linewidth, color='blue')

plt.xlabel('Current amplitude [µA]', fontsize=20)
plt.ylabel('Nb. of active neurons', fontsize=20)
plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40], fontsize=15)
plt.yticks(fontsize=15)
plt.title('Activity for increasing stim. amplitude', fontsize=20)
plt.tight_layout()



############################################################
############ DATA ANALYSIS DEPTH AND LOCAL #################
############################################################

nbs = []
centroids = []

for slice in list_depth:
    print(slice)
    for pattern in range(0,slice.nb_patterns):
        print('Pattern: ' + str(pattern))
        nb_active_neurons = slice.get_nb_active_neurons(pattern=pattern)[0]
        _, centroid, stdev, _, = slice.spatial_analysis(pattern=pattern, fitting_algorithm='gaussian', plot=False)
        nbs.append(nb_active_neurons)
        centroids.append(centroid)

plt.figure()
markersize=7
colors = ['blue', 'orange', 'blue', 'orange', 'blue', 'orange', 'blue', 'orange']

for centroid in range(0, len(centroids), 4):
    plt.plot([1, 2, 3, 4], [centroids[centroid], centroids[centroid+1], centroids[centroid+2], centroids[centroid+3]], color='black', alpha=0.5)

for centroid in range(0, len(centroids), 8):
    plt.plot(1, centroids[centroid], 'o', markersize=markersize, color=colors[0], label='20 µA')
for centroid in range(4, len(centroids), 8):
    plt.plot(1, centroids[centroid], 'o', markersize=markersize, color=colors[1], label='10 µA')

for centroid in range(1, len(centroids), 8):
    plt.plot(2, centroids[centroid], 'o', markersize=markersize, color=colors[2])
for centroid in range(5, len(centroids), 8):
    plt.plot(2, centroids[centroid], 'o', markersize=markersize, color=colors[3])

for centroid in range(2, len(centroids), 8):
    plt.plot(3, centroids[centroid], 'o', markersize=markersize, color=colors[4])
for centroid in range(6, len(centroids), 8):
    plt.plot(3, centroids[centroid], 'o', markersize=markersize, color=colors[5])

for centroid in range(3, len(centroids), 8):
    plt.plot(4, centroids[centroid], 'o', markersize=markersize, color=colors[6])
for centroid in range(7, len(centroids), 8):
    plt.plot(4, centroids[centroid], 'o', markersize=markersize, color=colors[7])

plt.ylim([-0.1, 1.1])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0\nElec. A', '30', '60', '90', '120', 'Elec. B\n150'], fontsize=15)
plt.ylabel('Location of centroid \nalong layer [µm]', fontsize=20)
plt.xticks([1,2,3,4], ['A to B', 'A to low B', 'Low A to B', 'Low A & B'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Centroid shift for depth stim.', fontsize=20)
plt.tight_layout()


plt.figure()
markersize=7
colors = ['blue', 'orange', 'blue', 'orange', 'blue', 'orange']

for centroid in range(0, len(nbs), 4):
    plt.plot([1, 2, 3], [nbs[centroid], (nbs[centroid+1] + nbs[centroid+2])/2, nbs[centroid+3]], color='black', alpha=0.5)

for centroid in range(0, len(nbs), 8):
    plt.plot(1, nbs[centroid], 'o', markersize=markersize, color=colors[0], label='20 µA')
for centroid in range(4, len(nbs), 8):
    plt.plot(1, nbs[centroid], 'o', markersize=markersize, color=colors[1], label='10 µA')

for centroid in range(1, len(nbs), 8):
    plt.plot(2, (nbs[centroid] + nbs[centroid+1])/2, 'o', markersize=markersize, color=colors[2])
for centroid in range(5, len(nbs), 8):
    plt.plot(2, (nbs[centroid] + nbs[centroid+1])/2, 'o', markersize=markersize, color=colors[3])

for centroid in range(3, len(nbs), 8):
    plt.plot(3, nbs[centroid], 'o', markersize=markersize, color=colors[4])
for centroid in range(7, len(nbs), 8):
    plt.plot(3, nbs[centroid], 'o', markersize=markersize, color=colors[5])

plt.yticks(fontsize=15)
plt.ylabel('Nb. of activated neurons', fontsize=20)
plt.xticks([1,2,3], ['A to B in-plane', '1 lower electrode', '2 lower electrodes'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Amount of activity for depth stim.', fontsize=20)
plt.tight_layout()




nbs = []
centroids = []
stdevs = []

for slice in list_local:
    print(slice)
    for pattern in range(0,slice.nb_patterns):
        print('Pattern: ' + str(pattern))
        nb_active_neurons = slice.get_nb_active_neurons(pattern=pattern)[0]
        _, centroid, stdev, _, = slice.spatial_analysis(pattern=pattern, fitting_algorithm='gaussian', plot=False)
        nbs.append(nb_active_neurons)
        centroids.append(centroid)
        stdevs.append(stdev)

plt.figure()
markersize=7
colors = ['blue', 'orange', 'blue', 'orange', 'blue', 'orange']

for centroid in range(0, len(centroids), 3):
    plt.plot([1, 2, 3], [centroids[centroid], centroids[centroid+1], centroids[centroid+2]], color='black', alpha=0.5)

for centroid in range(0, len(centroids), 6):
    plt.plot(1, centroids[centroid], 'o', markersize=markersize, color=colors[0], label='20 µA')
for centroid in range(3, len(centroids), 6):
    plt.plot(1, centroids[centroid], 'o', markersize=markersize, color=colors[1], label='10 µA')

for centroid in range(1, len(centroids), 6):
    plt.plot(2, centroids[centroid], 'o', markersize=markersize, color=colors[2])
for centroid in range(4, len(centroids), 6):
    plt.plot(2, centroids[centroid], 'o', markersize=markersize, color=colors[3])

for centroid in range(2, len(centroids), 6):
    plt.plot(3, centroids[centroid], 'o', markersize=markersize, color=colors[4])
for centroid in range(5, len(centroids), 6):
    plt.plot(3, centroids[centroid], 'o', markersize=markersize, color=colors[5])

plt.ylim([-0.1, 1.1])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0\nElec. A', '30', '60', '90', '120', 'Elec. B\n150'], fontsize=15)
plt.ylabel('Location of centroid \nalong layer [µm]', fontsize=20)
plt.xticks([1,2,3], ['Regular (A to B)', 'Hyper-local (A)', 'Hyper-local (B)'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Centroid shift for hyper-local stim.', fontsize=20)
plt.tight_layout()


plt.figure()
markersize=7
colors = ['blue', 'orange', 'blue', 'orange', 'blue', 'orange']

for centroid in range(0, len(nbs), 3):
    plt.plot([1, 2], [nbs[centroid], (nbs[centroid+1] + nbs[centroid+2])/2], color='black', alpha=0.5)

for centroid in range(0, len(nbs), 6):
    plt.plot(1, nbs[centroid], 'o', markersize=markersize, color=colors[0], label='20 µA')
for centroid in range(3, len(nbs), 6):
    plt.plot(1, nbs[centroid], 'o', markersize=markersize, color=colors[1], label='10 µA')

for centroid in range(1, len(nbs), 6):
    plt.plot(2, (nbs[centroid] + nbs[centroid+1])/2, 'o', markersize=markersize, color=colors[2])
for centroid in range(4, len(nbs), 6):
    plt.plot(2, (nbs[centroid] + nbs[centroid+1])/2, 'o', markersize=markersize, color=colors[3])


plt.yticks(fontsize=15)
plt.ylabel('Nb. of activated neurons', fontsize=20)
plt.xticks([1,2], ['Regular bipolar', 'Hyper-local'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Amount of activity for hyper-local stim.', fontsize=20)
plt.tight_layout()


plt.figure()
markersize=7
colors = ['blue', 'orange', 'blue', 'orange', 'blue', 'orange']

for centroid in range(0, len(nbs), 3):
    plt.plot([1, 2], [stdevs[centroid], (stdevs[centroid+1] + stdevs[centroid+2])/2], color='black', alpha=0.5)

for centroid in range(0, len(nbs), 6):
    plt.plot(1, stdevs[centroid], 'o', markersize=markersize, color=colors[0], label='20 µA')
for centroid in range(3, len(nbs), 6):
    plt.plot(1, stdevs[centroid], 'o', markersize=markersize, color=colors[1], label='10 µA')

for centroid in range(1, len(nbs), 6):
    plt.plot(2, (stdevs[centroid] + stdevs[centroid+1])/2, 'o', markersize=markersize, color=colors[2])
for centroid in range(4, len(nbs), 6):
    plt.plot(2, (stdevs[centroid] + stdevs[centroid+1])/2, 'o', markersize=markersize, color=colors[3])


plt.yticks(fontsize=15)
plt.ylabel('Std. dev. of neuronal spread', fontsize=20)
plt.xticks([1,2], ['Regular bipolar', 'Hyper-local'], fontsize=15)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles)) # To remove duplicate labels because of the for loop
plt.legend(by_label.values(), by_label.keys(), fontsize=12)

plt.title('Extent of activity for hyper-local stim.', fontsize=20)
plt.tight_layout()

plt.show()