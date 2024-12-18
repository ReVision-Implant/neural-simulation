import sys
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path)
from plotter_helper import Plotter
import matplotlib.pyplot as plt

exp = [2]
patterns = [7,8]
amplitudes=[10]
mice=[0,1]
n_rows=2
n_cols=2
row_param ="mice"
col_param ="patterns"
dir='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_2/output/plots/spatial/layer4_symm_ampl.png'

plot=Plotter(n_rows,n_cols)



plot.plot_all(exp, patterns, amplitudes, mice, row_param, col_param, depth=None, sigma=10, centroid=True, electrodes=True, ticks=True)
plot.legend()
titles = [['layer 4 symm, 10 µA', 'layer 4 asymm, 20 µA'], ['layer 4 symm, 20 µA', 'layer 4 asymm, 20µA']]
plot.set_titles(titles)

#plt.savefig(dir)
plt.show()
