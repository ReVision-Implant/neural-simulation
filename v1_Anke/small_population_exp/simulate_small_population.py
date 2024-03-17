#This script builds a small population of 100 neurons. The goal is to test different types of axons and parameters
#Option: make axons longer (adjust code run file)

from bmtk.utils.sim_setup import build_env_bionet
from bmtk.simulator import bionet


build_env_bionet(
    base_dir='simulation_axon_10_diam_1',
    #base_dir='simulation_long_4',
    config_file='config.json',
    network_dir='network_mask_dense',
    tstop=100.0, dt=0.025,
    #report_vars=['v'], # Record membrane potential
    #current_clamp={  # Creates a step current from 500.0 ms to 1500.0 ms  
     #   'amp': 0.120,
    #    'delay': 500.0,
    #    'duration': 1000.0
    #},
    #include_examples=True,    # Copies components files
    #compile_mechanisms=True #will try to compile NEURON mechanisms -> often fails and must be done manually; left it in so it will give an error and remember me to do it manually

)