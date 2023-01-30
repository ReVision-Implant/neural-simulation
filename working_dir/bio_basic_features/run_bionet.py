"""Simulates an example network of 450 cell receiving two kinds of exernal input as defined in the configuration file"""
import sys
import matplotlib.pyplot as plt
import h5py
import numpy as np
import json
import os 

from bmtk.simulator import bionet
from bmtk.analyzer.compartment import plot_traces
from biophys_components.voltage_waveform import CreateVoltageWaveform

def show_cell_var(conf, var_name):
    plot_traces(config_file=conf, report_name=var_name)


def run(config_file):
    
    with open (config_file, "r") as f:
        parameters = json.load(f)

    time_step = parameters["run"]["time_step"]
    currentamplitude = parameters["run"]["current_amplitude"]

    print("time_step =", time_step)
    print("current_amplitude =", currentamplitude)

    #creates the file xstim
    CreateVoltageWaveform(current_amplitude=currentamplitude, timestep=time_step, plotting=True) 

    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()

    show_cell_var(config_file, 'membrane_potential')


# CreateVoltageWaveform(current_amplitude=20, timestep=10, plotting=True)



if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        # Make sure to run only one at a time
        run('config_iclamp.json')  # Current clamp stimulation
        # run('config_xstim.json')  # Extracellular electrode stimulation
        # run('config_spikes_input.json')  # Synaptic stimulation with external virtual cells
