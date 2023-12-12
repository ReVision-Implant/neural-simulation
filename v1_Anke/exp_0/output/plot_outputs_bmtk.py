config_file_test='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_0/config/config_test_1.json'
spike_file_test='/scratch/leuven/356/vsc35693/neural-simulation/v1_Anke/exp_0/output/pattern_0/amplitude_10/test_nomask_1/spikes.h5'

from bmtk.analyzer.spike_trains import to_dataframe

results_df=to_dataframe(config_file_test,spike_file_test,population='v1')
print('Number of Spikes mouse 1 no mask pattern 0: {}'.format(len(results_df)))
results_df.head()

#from bmtk.analyzer.spike_trains import plot_raster, plot_rates_boxplot

#plot_raster(config_file_test, group_by='pop_name')
#plot_rates_boxplot(config_file=config_file_test, group_by='pop_name')