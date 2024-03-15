from bmtk.analyzer.spike_trains import plot_raster

number_stubs='10'
diam='5'
conductance='X'
title_plot='simulation 10  axon stubs basic'

plot_raster(spikes_file= 'simulation_long_10/output/spikes.csv', title=title_plot)
#plot_raster(config_file = 'simulation_long_'+number_stubs+'_diam_'+diam+'/config.json', spikes_file='simulation_long_'+number_stubs+'_diam_'+diam+'/output/spikes.csv', title=title_plot)
