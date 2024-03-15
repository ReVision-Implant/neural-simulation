from bmtk.analyzer.spike_trains import plot_raster

number_stubs='10'
diam='5'
title_plot='simulation no axons'

plot_raster(config_file = 'simulation/config.json', spikes_file= 'simulation/output/spikes.csv', title=title_plot)
#plot_raster(config_file = 'simulation_long_'+number_stubs+'_diam_'+diam+'/config.json', spikes_file='simulation_long_'+number_stubs+'_diam_'+diam+'/output/spikes.csv', title=title_plot)