import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import norm
from file_helper import format_params, get_dirs
from spikes_helper import get_grid, get_spikes, get_centroid_cov

class SubPlotter():
    """The SubPlotter class can be used together with the Plotter class to modify individual subplots. 
    In most cases, it is probably easier to use Plotter's plot_all(), which automatically calls this class.
    """

    def __init__(self, ax, exp, pattern, amplitude, mice, depth=None, sigma=10, centroid=True, electrodes=True, ticks=True):
        """Setting attributes and calling main functions

        :param ax: (maplotlib Axes object) Axes of the subplot.
        :param exp: (int/str) Experiment name.
        :param patterns: (int/str or list thereof) Pattern names. 
        :param amplitudes: (int/str or list thereof) Amplitudes in uA. 
        :param mice: (int/str or list thereof) Mouse names. 
        :param depth: (tuple or list with 2 elements) If not None, only include spikes of neurons that are in the layer between depth[0] and depth[1]. Defaults to None
        :param sigma: (int) Standard deviation of the Gaussian smoother. Defaults to 10.
        :param centroid: (bool) If True, plot centroid. Defaults to True.
        :param electrodes: (bool) If True, plot electrodes. Defaults to True.
        :param ticks: (bool) If True, set ticks to default value. Defaults to True
        """   
        ### Setting attributes
        self.ax = ax
        exp, pattern, amplitude, mice = format_params(exp, pattern, amplitude, mice)
        self.exp = exp
        self.pattern = pattern[0]
        self.amplitude = amplitude[0]
        self.mice = mice
        self.depth = depth
        self.sigma = sigma

        ### Main functions
        self.set_dirs()
        self.get_image()
        self.plot_image()
        if centroid: self.plot_centroid()
        if electrodes: self.plot_electrodes()
        if ticks: self.set_ticks()

    def set_dirs(self):
        """Set directories as class attributes""" 
        kwargs = get_dirs(self.exp, self.pattern, self.amplitude, self.mice)
        for k,v in kwargs.items():
            setattr(self, k, v)
        
        return

    def get_image(self):
        """Create an intensity map of neural activity based on the number of spikes.
        Spikes are binned into 1 um^2 pixels and smoothed with a Gaussian filter.
        The area of the square plot outside of the neural column is set to NaN.
        """
        ### Initalisation
        self.size = int(self.radius*2)
        self.grid = np.zeros((self.size,self.size))

        ### Get the node positions and corresponding spike count of each neuron
        self.node_pos, self.n_spikes = get_spikes(**self.__dict__)

        ### Get grid based on spikes
        self.grid = get_grid(**self.__dict__)

        ### Gaussian smoothing of the grid
        self.image = gaussian_filter(self.grid, self.sigma, truncate=4) * self.sigma / len(self.nodes_dirs)

        ### Set colour outside circle to NaN
        for i in range(len(self.image)):
            for j in range(len(self.image)):
                if (i-self.size/2+.5)**2 + (j-self.size/2+.5)**2 >= (self.size/2)**2:
                    self.image[i,j] = float('NaN')

        return
    

    def plot_image(self):
        """Plot the intensity map of neural activity created by get_image().
        The intensities are converted to an image with a green color map.
        Inactive background is set to black and the NaN values are set to white.
        """        
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
    

    def plot_electrodes(self, **kwargs):
        """Plot the electrode locations as specified in pattern.csv.
        Electrodes with positive and negative amplitudes are plotted in a different colour.

        :param **kwargs: Passed to ax.scatter().
        """
        
        ### Get electrode names from pattern.csv
        pattern_dir = self.pattern_dirs[0]
        pattern = pd.read_csv(pattern_dir, sep='\s+')
        electrode_names = pattern['electrode'].to_numpy()

        ### Get electrode positions from electrodes.csv
        electrodes_dir = self.electrodes_dirs[0]
        electrodes = pd.read_csv(electrodes_dir, sep='\s+')
        electrodes = electrodes.loc[electrodes['electrode'].isin(electrode_names)]
        self.electrode_pos = electrodes[['pos_x','pos_y','pos_z']].to_numpy()

        ### Get amplitudes from pattern.csv
        amplitudes = pattern['amplitude'].to_numpy()
        elec_pos = self.electrode_pos[amplitudes>0, :]
        elec_neg = self.electrode_pos[amplitudes<0, :]

        ### Plotting
        self.ax.scatter(elec_pos[:,0], elec_pos[:,2], marker='X', s=200, c='red', edgecolors='black', label='stimulating\nelectrode', alpha=0.7, **kwargs)
        self.ax.scatter(elec_neg[:,0], elec_neg[:,2], marker='X', s=200, c='white', edgecolors='black', label='return\nelectrode', alpha=0.7, **kwargs)
    
        return


    def plot_centroid(self, **kwargs):
        """Get the centroid location and plot it. 
        
        :param **kwargs: Passed to ax.scatter().
        """        
        ### Get centroid
        self.centroid, self.cov = get_centroid_cov(self.node_pos, self.n_spikes)

        ### Plotting
        self.ax.scatter(self.centroid[0], self.centroid[2], marker='X', s=200, c='lime', edgecolors='black', label='centroid', alpha=0.7, **kwargs)

        return    
    

    def show_ratio(self):
        """Show the ratio between the centroid-to-electrode distances.
        Only makes sense in cases where two electrodes stimulate with opposite polarity.
        """    
        ### Show centroid-to-electrode distance ratio |centroid->electrode0|/(|centroid->electrode0| + |centroid->electrode1|)
        if np.size(self.electrode_pos,0) > 1:
            c = self.centroid
            e = self.electrode_pos
            ratio = np.linalg.norm(c[[0,2]] - e[0,[0,2]]) / (np.linalg.norm(c[[0,2]] - e[1,[0,2]]) + np.linalg.norm(c[[0,2]] - e[0,[0,2]]))
        self.ax.set_xlabel('centroid ratio = ' + str(np.round(ratio,2)))

        return
    
    def set_ticks(self, xticks=[-182,0,182], yticks=[-182,0,182], **kwargs):
        """Set axes ticks. Defaults to the 3x3 electrode grid locations

        :param xticks: Passed to ax.set_xticks(). Defaults to [-182,0,182].
        :param yticks: Passed to ax.set_yticks(). Defaults to [-182,0,182].
        :param **kwargs: Passed to ax.tick_params().
        """
        self.ax.set_xticks(xticks)
        self.ax.set_yticks(yticks)
        self.ax.tick_params(direction='inout', bottom=True, top=True, left=True, right=True, width=3, length=7, **kwargs)

        return
    

class Plotter():
    """Main plotter class that creates a figure consisting of a number of subplots.
    Each subplot is generated and modified with the SubPlotter class.
    However, plot_all() and set_all() can be used to generate or modify all subplots at once.
    """

    def __init__(self, n_rows=1, n_cols=1):
        """Set attributes and create a fig and axs object for the subplots

        :param n_rows: (int) Number of rows of the subplot grid. Defaults to 1.
        :param n_cols: (int) Number of columns of the subplot grid. Defaults to 1.
        """   
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.fig, self.axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(5*n_cols,5*n_rows))

    def plot_all(self, exp, patterns, amplitudes, mice, row_param=None, col_param=None, depth=None, sigma=10, centroid=True, electrodes=True, ticks=True):
        """Iteratively populate all subplots. For more control over each specific subplot, use SubPlotter explicitly.
        All parameters (exp, patterns, amplitudes, mice) can be varied between subplot rows and columns.
        When there is only one row/column, row_param/col_param should be None.
        If row_param/col_param is not 'mice', the images of all specified mice are combined in all subplots.
        E.g.
        plot = Plotter(n_rows=3, n_cols=2)
        plot.plot_all(exp=[1,2], patterns=[0,1,2], amplitudes=20, mice=[0,1,2], row_param='patterns', col_param='exp')
        will create a 3x2 grid of plots where each plot shows the results of (exp[col], pattern[row], amplitude, mice[:])
        E.g.
        plot = Plotter(n_rows=1, n_cols=3)
        plot.plot_all(exp=1, patterns=1, amplitudes=10, mice=[0,1,2], row_param=None, col_param='mice')
        will create a 1x3 grid of plots where each plot shows the results of (exp, pattern, amplitude, mice[col])

        :param exp: (int/str) Experiment name.
        :param patterns: (int/str or list thereof) Pattern names. 
        :param amplitudes: (int/str or list thereof) Amplitudes in uA. 
        :param mice: (int/str or list thereof) Mouse names.
        :param row_param: (str or None) If not None, specifies the parameter that is varied between subplot rows. Else, n_rows should be 1. Defaults to None.
        :param col_param: (str or None) Idem but for subplot columns.
        :param depth: (tuple or list) If not None, only include spikes of neurons that are in the layer between depth[0] and depth[1]. Defaults to None
        :param sigma: (int) Standard deviation of the Gaussian smoother. Defaults to 10.
        :param centroid: (bool) If True, plot centroid. Defaults to True.
        :param electrodes: (bool) If True, plot electrodes. Defaults to True.
        :param ticks: (bool) If True, set ticks to default value. Defaults to True
        """        
        exp, patterns, amplitudes, mice = format_params(exp, patterns, amplitudes, mice)    
        self.exp = exp
        self.patterns = patterns
        self.amplitudes = amplitudes
        self.mice = mice
        self.depth = depth
        self.sigma = sigma

        rows = getattr(self, row_param) if row_param is not None else [None]
        cols = getattr(self, col_param) if col_param is not None else [None]

        self.subplots = np.empty((self.n_rows,self.n_cols), dtype=object)

        assert len(rows) == self.n_rows and len(cols) == self.n_cols

        for row, row_value in enumerate(rows):
            for col, col_value in enumerate(cols):
                if row_value is not None: setattr(self, row_param, row_value)
                if col_value is not None: setattr(self, col_param, col_value)
                if len(rows) == 1 and len(cols) == 1:
                    ax = self.axs
                elif len(rows) == 1 or len(cols) == 1:
                    ax = self.axs[row+col]
                else:
                    ax = self.axs[row, col]
                self.subplots[row,col] = SubPlotter(ax, self.exp, self.patterns, self.amplitudes, self.mice, self.depth, self.sigma, centroid=centroid, electrodes=electrodes, ticks=ticks)
    
        return
    

    def set_all(self, func, **kwargs):
        """Calls a function from the SubPlotter class for every SubPlotter object. 
        E.g. 
        :param func: (str) Name of the SubPlotter class function.
        :param **kwargs: Passed to SubPlotter.func(**kwargs)
        """
        for i in range(len(self.subplots.flatten())):
            subplot = self.subplots.flatten()[i]
            getattr(subplot, func)(**kwargs)
                
        return
    
    def legend(self, **kwargs):
        """Sets the legend in the figure.
        First collects the handles and labels of all elements in any subplot.
        Then removes elements with duplicate label.
        fig = Plotter(1,2)
        fig.plot_all(0,0,30,[1,2],'patterns','mice', False, False, False)
        fig.set_all(func='set_ticks')
        fig.set_all(func='plot_centroid')
        fig.set_all(func='plot_electrodes')
        :param **kwargs: Passed to self.fig.legend(**kwargs)
        """        
        handles, labels = zip(*[ax.get_legend_handles_labels() for ax in self.axs])
        handles = [y for x in handles for y in x]
        labels = [y for x in labels for y in x]
        i = 0
        while i < len(labels):
            if labels[i] in labels[:i]:
                labels[i:i+1] = []
                handles[i:i+1] = []
                i -= 1
            i += 1
        self.fig.legend(handles=handles, labels=labels, **kwargs)


if __name__ == "__main__":

    plt.rcParams.update({'font.size':30})
    plt.rcParams.update({'legend.fontsize':30})
    plt.rcParams.update({'axes.labelsize':30})
    plt.rcParams['font.family'] = 'Times New Roman'
    mfont = {'fontname':'Times New Roman'}

    fig = Plotter(1,2)
    fig.plot_all(0,0,30,[1,2],'patterns','mice')
    fig.set_all(func='set_ticks')
    fig.set_all(func='plot_centroid')
    fig.set_all(func='plot_electrodes')

    # fig.legend()

    # fig = Plotter(2,2)
    # test = SubPlotter(fig.axs[0,0], 0, 0, 30, [1,2])
    # test = SubPlotter(fig.axs[0,1], 0, 0, 30, 2)
    # test = SubPlotter(fig.axs[1,0], 0, 0, 30, 1)
    # test = SubPlotter(fig.axs[1,1], 0, 0, 30, [1,2])
    
    plt.show()