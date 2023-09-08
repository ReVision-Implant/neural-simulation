import pandas as pd
from hdf5 import HDF5
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import norm
from file_helper import format_params, get_dirs
from spikes_helper import get_grid, get_spikes, get_centroid_cov

class SubPlotter():

    def __init__(self, ax, exp, pattern, amplitude, mice, depth=None, sigma=10):

        self.ax = ax

        exp, pattern, amplitude, mice = format_params(exp, pattern, amplitude, mice)
        self.exp = exp
        self.pattern = pattern[0]
        self.amplitude = amplitude[0]
        self.mice = mice

        self.depth = depth
        self.sigma = sigma
        self.set_dirs()
        self.node_pos, self.n_spikes = get_spikes(**self.__dict__)
        self.get_image()
        self.plot_image()
        self.plot_centroid()
        self.plot_electrodes()
        self.set_ticks()

    def set_dirs(self):
        
        kwargs = get_dirs(self.exp, self.pattern, self.amplitude, self.mice)
        for k,v in kwargs.items():
            setattr(self, k, v)
        
        return

    # def get_spikes(self):
        
    #     self.node_pos, self.n_spikes = get_spikes(**self.kwargs)

    #     return

    def get_image(self):

        ### Initalisation
        self.size = int(self.radius*2)
        self.grid = np.zeros((self.size,self.size))

        ### Get spikes from each simulation and add them in a single grid array
        self.grid = get_grid(**self.__dict__)

        ### Gaussian smoothing of the image
        self.image = gaussian_filter(self.grid, self.sigma, truncate=4) * self.sigma / len(self.nodes_dirs)

        ### Set colour outside circle to NaN
        for i in range(len(self.image)):
            for j in range(len(self.image)):
                if (i-self.size/2+.5)**2 + (j-self.size/2+.5)**2 >= (self.size/2)**2:
                    self.image[i,j] = float('NaN')

        return
    

    def plot_image(self):
        
        ### Generate green color map with black background that also maps masked valued to white 
        cdict = {
            'red': (
                (0.0,  0.0, 0.0),
                (1.0,  0.0, 0.0),
            ),
            'green': (
                (0.0,  0.0, 0.0),
                (0.1,  0.3, 0.3),
                (0.2,  0.45, 0.45),
                (0.4,  0.66, 0.66),
                (0.6,  0.8, 0.8),   
                (0.8,  0.9, 0.9),
                (1.0,  1.0, 1.0),
            ),
            'blue': (
                (0.0,  0.0, 0.0),
                (1.0,  0.0, 0.0),
            )
        }
        cmap = matplotlib.colors.LinearSegmentedColormap('green', cdict)
        cmap.set_bad(color='white')

        ### Create image mask from NaN values (outside circle) and plot the resulting masked image
        masked_array = np.ma.masked_where(np.isnan(self.image), self.image)
        im = self.ax.imshow(masked_array, cmap=cmap, extent=[-self.radius,self.radius,-self.radius,self.radius], origin = 'lower', vmin=0, vmax=1)

        return

    def plot_electrodes(self):

        ### Plotting options
        s=200
        alpha = 0.7
        marker = 'X'
        c_pos = 'red'
        c_neg = 'white'
        edgecolors = 'black'
        label_pos = 'stimulating\nelectrode'
        label_neg = 'return\nelectrode'

        ### Extract electrode positions and polarity from pattern.csv
        pattern_dir = self.pattern_dirs[0]
        pattern = pd.read_csv(pattern_dir, sep=' ')
        self.electrode_pos = pattern[['pos_x', 'pos_y', 'pos_z']].to_numpy()
        polarity = pattern['amplitude'].to_numpy()
        elec_pos = self.electrode_pos[polarity>0, :]
        elec_neg = self.electrode_pos[polarity<0, :]

        ### Plotting
        self.ax.scatter(elec_pos[:,0], elec_pos[:,2], marker=marker, s=s, c=c_pos, edgecolors=edgecolors, label=label_pos, alpha=alpha)
        self.ax.scatter(elec_neg[:,0], elec_neg[:,2], marker=marker, s=s, c=c_neg, edgecolors=edgecolors, label=label_neg, alpha=alpha)
    
        return


    def plot_centroid(self):

        ### Plotting options
        s = 200
        marker = 'X'
        c = 'lime'
        edgecolors = 'black'
        label = 'centroid'
        alpha = 0.7

        ### Get centroid
        self.centroid, self.cov = get_centroid_cov(self.node_pos, self.n_spikes)

        ### Plotting
        self.ax.scatter(self.centroid[0], self.centroid[2], marker=marker, s=s, c=c, edgecolors=edgecolors, label=label, alpha=alpha)

        return    

    def show_ratio(self):

        ### Show centroid-to-electrode distance ratio |centroid->electrode0|/(|centroid->electrode0| + |centroid->electrode1|)
        if np.size(self.electrode_pos,0) > 1:
            c = self.centroid
            e = self.electrode_pos
            ratio = np.linalg.norm(c[[0,2]] - e[0,[0,2]]) / (np.linalg.norm(c[[0,2]] - e[1,[0,2]]) + np.linalg.norm(c[[0,2]] - e[0,[0,2]]))
        self.ax.set_xlabel('centroid ratio = ' + str(np.round(ratio,2)))

        return
    
    def set_ticks(self):
        
        if self.radius == 200:
            self.ax.set_xticks([-91,91])
            self.ax.set_yticks([-91,91])
        elif self.radius == 400:
            self.ax.set_xticks([-182,0,182])
            self.ax.set_yticks([-182,0,182])
        self.ax.tick_params(direction='inout', bottom=True, top=True, left=True, right=True, width=3, length=7)

        return
    

class Plotter():

    def __init__(self, n_rows, n_cols):
        
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.fig, self.axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(4*n_cols,4*n_rows))

    def plot_all(self, exp, patterns, amplitudes, mice, row_param=None, col_param=None):
        
        exp, patterns, amplitudes, mice = format_params(exp, patterns, amplitudes, mice)    
        self.exp = exp
        self.patterns = patterns
        self.amplitudes = amplitudes
        self.mice = mice

        rows = getattr(self, row_param)
        cols = getattr(self, col_param)

        assert len(rows) == self.n_rows and len(cols) == self.n_cols

        for row, row_value in enumerate(rows):
            for col, col_value in enumerate(cols):
                setattr(self, row_param, row_value)
                setattr(self, col_param, col_value)
                if len(rows) == 1 and len(cols) == 1:
                    SubPlotter(self.axs, self.exp, self.patterns, self.amplitudes, self.mice)
                elif len(rows) == 1 or len(cols) == 1:
                    SubPlotter(self.axs[row+col], self.exp, self.patterns, self.amplitudes, self.mice)
                else:
                    SubPlotter(self.axs[row,col], self.exp, self.patterns, self.amplitudes, self.mice)
    
    def set_all_ticks(self):

        for row in self.n_rows:
            for col in self.n_cols:
                self.axs[row,col].set_xticks([-182,0,182])
                self.axs[row,col].set_yticks([-182,0,182])
                self.axs[row,col].tick_params(direction='inout', bottom=True, top=True, left=True, right=True, width=3, length=7)

        return
    


if __name__ == "__main__":

    plotter_test = Plotter(1,2)
    plotter_test.plot_all(0,0,30,[1,2],'patterns','mice')

    subplotter_test = Plotter(2,2)
    test = SubPlotter(subplotter_test.axs[0,0], 0, 0, 30, [1,2])
    test = SubPlotter(subplotter_test.axs[0,1], 0, 0, 30, 2)
    test = SubPlotter(subplotter_test.axs[1,0], 0, 0, 30, 1)
    test = SubPlotter(subplotter_test.axs[1,1], 0, 0, 30, [1,2])
    plt.show()

