import sys
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path)
from plotter_helper import Plotter

exp = [2]
patterns = [0,4,5,6,7,8]
amplitudes=[10,20]
mice=[0,1]
n_rows=6
n_cols=2
row_param ="patterns"
col_param ="mice"

plot=Plotter(n_rows,n_cols)

plot.plot_all(exp, patterns, amplitudes, mice, row_param, col_param, depth=None, sigma=10, centroid=True, electrodes=True, ticks=True)