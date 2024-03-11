#This script builds a small population of 100 neurons. The goal is to test different types of axons and parameters
#Option 1: changing the direction of the axon with aibs_perisomatic_directed (change in build file first!)
#Option 2: make axons longer (adjust code run file)

from bmtk.utils.sim_setup import build_env_bionet
from bmtk.simulator import bionet


build_env_bionet(
    X#base_dir='simulation',
    #base_dir='simulation_long_axons',
    #base_dir='simulation_directed_axons', #use directed axons
    config_file='config.json',
    network_dir='network',
    #network_dir='network_directed_axons', #use the directed axons
    tstop=2000.0, dt=0.1,
    report_vars=['v'], # Record membrane potential
    current_clamp={  # Creates a step current from 500.0 ms to 1500.0 ms  
        'amp': 0.120,
        'delay': 500.0,
        'duration': 1000.0
    },
    include_examples=True,    # Copies components files
    compile_mechanisms=True #will try to compile NEURON mechanisms -> often fails and must be done manually; left it in so it will give an error and remember me to do it manually

)