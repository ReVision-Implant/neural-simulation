import sys;
module_path='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/toolbox';
sys.path.append(module_path);

from plotter_helper import Plotter;

plot=Plotter(n_rows=1,n_cols=2);
plot.plot_all(exp=0, patterns=[0], amplitudes=10, mice=[0,1],row_param=None,col_param='mice');