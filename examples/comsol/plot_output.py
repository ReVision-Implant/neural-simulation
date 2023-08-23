import matplotlib.pyplot as plt
import sys
sys.path.append('../../..')
sys.path.append('../..')
sys.path.append('..')
sys.path.append('../bio_components')
from bmtk.analyzer.compartment import plot_traces
from bmtk.analyzer.spike_trains import plot_raster, plot_rates_boxplot

plot_raster(config_file='config_comsol_1c.json')
plot_traces(config_file='config_comsol_1c.json', report_name='membrane_potential')


