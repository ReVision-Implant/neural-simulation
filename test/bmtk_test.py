import bmtk
import bmtk.analyzer.compartment as comp

print(comp.__file__)
comp.plot_traces(config_file='$VSC_SRCATCH/neural-simulation/examples/bio_comsol/config.comsol_stat.json', save_as='$VSC_SCRATCH/neural-simulation/test/save')
print('OK')