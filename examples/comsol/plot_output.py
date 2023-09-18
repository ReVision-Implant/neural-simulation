import matplotlib.pyplot as plt

from bmtk.analyzer.compartment import plot_traces
from bmtk.analyzer.spike_trains import plot_raster

plot_raster(config_file='config.comsol_tdep.json')
plot_traces(config_file='config.comsol_tdep.json', report_name='membrane_potential')

plot_raster(config_file='config.comsol_stat.json')
plot_traces(config_file='config.comsol_stat.json', report_name='membrane_potential')

plot_raster(config_file='config.comsol_stat2.json')
plot_traces(config_file='config.comsol_stat2.json', report_name='membrane_potential')

plt.show()
