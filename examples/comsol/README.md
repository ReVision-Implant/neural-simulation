# COMSOL

A small network of 100 multi-compartment biophysically detailed cells placed in a column with a 100 $\mu m$ radius and 100 $\mu m$ height . 
The network receives input in the form of extracellular potentials defined by a .txt file that is exported from COMSOL. 

Uses the BioNet simulator (requires NEURON).

## Running:

```
$ python run_bionet.py config_comsol_*.json
```

All three configurations simulate the same stimulation configuration, where two probes, with one electrode each, are inserted into a cylindrical piece of tissue. Biphasic symmetric pulses are sent from one electrode to the other, i.e. they have opposite polarity. 

- config.comsol_tdep.json uses a time-dependent COMSOL simulation.
- config.comsol_stat.json uses a single stationary COMSOL simulation. The time dependency is described by the waveform in *examples/bio_components/stimulations/waveform.csv*.
- config.comsol_stat2.json uses two COMSOL simulations, in which either of the electrodes is active. The output of both simulations are superimposed to get the actual extracellular potentials. The time dependency is described by the waveform in examples/bio_components/stimulations/waveform.csv.

The output files have already been generated in the *outputs* directory containing log, spike trains and recorded cell variables. Running the simulations again will overwrite the existing files.

## The Network:
The network files have already been built and stored as SONATA files in the *network/* directory. The bmtk Builder
script used to create the network files is *build_network.py*. To adjust the parameters and/or topology of the network
change this file and run:
```
$ python build_network.py
```
This will overwrite the existing files in the network directory. Note that there is some randomness in how the network
is built, so expect (slightly) different simulation results everytime the network is rebuilt.

## Simulation Parameters
Parameters to run the simulation, including run-time, inputs, recorded variables, and networks are stored in config_comsol.json and config_network.json and can modified with a text editor.

The COMSOL file path can be specified in config_comsol.json

```json
  "inputs": {
    "Extracellular_Stim": {
      "input_type": "lfp",
      "node_set": "all",
      "module": "comsol",
      "comsol_file": "$STIM_DIR/comsol.txt",
      "waveform": "$STIM_DIR/waveform.csv",
      "amplitude": 10
    }
  }
```
where *waveform* is optional. If it is specified, the COMSOL file should contain the output of a stationary study. If not, the COMSOL file should contain the output of a time-dependent study.

It is also possible to specifiy multiple comsol_files and waveforms, in which case these arguments should be passed as a list, i.e. 

```json
  "inputs": {
    "Extracellular_Stim": {
      "input_type": "lfp",
      "node_set": "all",
      "module": "comsol",
      "comsol_file": ["$STIM_DIR/comsol1.txt", "$STIM_DIR/comsol2.csv"],
      "waveform": ["$STIM_DIR/waveform1.csv", "$STIM_DIR/waveform2.csv"],
      "amplitude": 10
    }
  }
```

## Plotting

*plot_output.py* will generate figures showing the results. Some of them are already saved in the *figures* directory.
