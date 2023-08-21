[Back to ToC](/docs/manual/README.md)

# Simulating the neural response with BMTK

In order to simulate the neural activation of the tissue in response to the extracellular potentials, we need two things:

- The extracellular potentials which were [calculated in COMSOL](../comsol/solution.md).
- A computational model of the neural tissue in BMTK.

Then, we can use the [comsol](/examples/comsol/README.md) module in BMTK's BioNet to import the COMSOL output and simulate the behaviour of the neural network in response to the imposed extracellular potentials.

## Building a network

Thanks to the work of [Billeh et al.](https://doi.org/10.1016/j.neuron.2020.01.040), a model of the mouse V1 already exists. However, it is possible to adapt certain parameters (e.g. the size or shape of the patch of tissue) to fit the specific needs of an experiment. Changing neuron models or synaptic connections is also possible, but obviously it is a lot more complex and might disrupt the proper functioning of the model.

In the context of ReVision's computational experiments, we will create multiple network models &ndash;each with a different random seed&ndash; in order to perform the same experiment multiple times on different virtual 'animals'. While the neuron types are constrained to a certain layer of V1 and the general connection rules between and within layers are set, the random seed determines the positions of the neurons within their respective layer, as well as the exact synaptic connections that are made between individual neurons.  

Building a network uses several components/scripts that are found in the v1 folder. 
- v1/build_files/ 
- v1/components/
- v1/build_network.py - Simple model tweaks can probably be made here.


Calling build_network.py in the terminal or with a bash/job script (with a few optional arguments) will actually build the network including any possible changes you might have made. 
```
$ python build_network.py -o [output] --fraction [fraction] --rng-seed [seed]
```
- [output] - Output directory. Defaults to networks_rebuilt/network.
- [fraction] - Fraction of the total neurons the model should include. Defaults to 1.
- [seed] - Random seed to build the network, changing this will allow the creation of different virtual 'animals'. Defaults to 100.

For parallel computing, add ```mpirun -np [np]``` before the command above.
- [np] - number of cores to use.

Calling build_network.py multiple times with different random seeds (make sure to also set a different output directory!) will instantiate several networks that represent the different virtual 'animals'.

An example of a bash script calling this function can be found [here](/v1/build.sh). 

## (Only for stationary comsol studies) generating waveform.csv

Skip this section if you use a time-dependent comsol study.

Generating the waveform.csv file is done with [/v1/components/voltage_waveform.py](/v1/components/voltage_waveform.py). 

## Running a simulation

Once you have built one or several networks, you can run simulations with the previously built network(s). This also requires extracellular potentials that were [calculated in COMSOL](../comsol/solution.md). Depending on the stimulation parameters, the COMSOL output should be either stationary or time-dependent.

Configuring the comsol input for BMTK in the config.json file, will look something like this

```json
    "Extracellular_Stim": {
        "input_type": "lfp",
        "node_set": "all",
        "module": "comsol",
        "comsol_file": "$STIM_DIR/exp3/02-.txt",
        "waveform": "$STIM_DIR/waveform.csv",
        "amplitude": 10,
        "ip_method": "L"
    }
```
- input_type - Always "lfp".
- node_set - Used to filter which cells receive the input, but here it probably does not make sense to use anything besides "all".
- module - Always "comsol".
- comsol_file - /path/to/comsol.txt from stationary or time-dependent study
- waveform - /path/to/waveform.csv. Only supply this with stationary comsol.txt, remove the entire line for time-dependent comsol.txt.
- amplitude - Multiplication factor for waveform. If the amplitudes in waveform.csv are normalised to [0,1], this can be used to set the current amplitude. Defaults to 1. 
- ip_method - "NN" for nearest neighbour, "L" for linear interpolation (only for stationary comsol.txt). Defaults to "NN".
