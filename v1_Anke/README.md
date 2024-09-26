# Biophysically detailed mouse v1 model (with active axons).

This is the main folder and contains files you can use to run simulations with the biophysically detailed V1 model with thalamacortical (LGN) and background (BKG) inputs. The folder (and this readme) is based on the simulations/v1_biophysical folder from the [V1 mouse model](https://www.dropbox.com/sh/w5u31m3hq6u2x5m/AACpYpeWnm6s_qJDpmgrYgP7a?dl=0) of the Allen Institute, and includes the additions/changes to use it for the simulation of electrical stimulation.The main difference with the v1 folder is that this folder contains files to simulate models including active axons stubs. 

 It is designed to run using the Allen Institute [Brain Modeling Toolkit](https://github.com/AllenInstitute/bmtk); but the network, input and config files are in the [SONATA data format](https://github.com/AllenInstitute/sonata) for use with other tools that support SONATA.


## Requirements
- Python 3.6+
- NEURON 7.4+
- BMTK 0.0.8+

At the moment it also requires access to a HPC cluster, as the network is too large to instantiate in even the most powerful of desktop machines (by 2019 standards).
A cluster with 100 Xeon cores will typically take a few hours to run through the entire 3.0 second simulation. 


## File Structure

<!---
* network_precomputed_syns/ - Similar to above but with the synaptic locations and synaptic weights fully specified. In the network/*edges_types.csv the columns weight_function,
     target_sections and distance_range are used by bmtk to calculate where and how to place synapses between every source/target pair. In these files the synaptic locations
     and weights are stored in the edges.h5 files, more complient with the SONATA format.
-->

### 1. Build_files
 - scripts and properties for (re)building the V1 network
 - in contrast to the original build file in bmtk, the build file in v1_Anke adjusts the node positions of neurons making sure they cannot coincide with the location of the probes
### 2. Components
- parameters, morphology, model and mechanisms files used to instantiate individual cells and synapses, as well as python files containing helper functions.
- under mechanisms/modfiles you will find a modfile called mammalian_spike_Anke.mod. This modfile was added to the orginal mechanisms of the V1 model and contains the equations defining the ionchannels in the active axon models. 
- New modfiles can be added (and compiled) under mechanisms/modfiles
- the subfolder stimulation contains the comsolfiles,(with the calculated extracellular potentials), mask file, location of the electrodes, definition of the patterns and applied waveforms.

### 3. Exp folders
- Each folder has minimum two subfolders: the config folder and the output folder. Experiments 2 and 4  contain  extra subfolders with plots of results.
  - config folder: contains the config files
  - output folder: sonata output files with spike rates
  - plot_layers: contains plots of the kernel density estimates for the neurons projected onto the cortical layers (parallel with the surface of the mouse brain)
  - plot_columns: kernel density estimate plots for the neurons projected perpendicular to the layers (a coronal view on the cortex)

- The numbers of the experiments refer to:
  - exp_2: passive axon stubs as in original
  - exp_3: passive axon stubs + no virtual connections (synaptic weights of the neurons were set to zero)
  - exp_4: active axon stubs
    - This folder contains two extra plot folders: The spearman plots plot the spike rates of different experiments against eachother. The compare_exp2 plots pattern 0 for the model with and without axons next to eachother.
  - exp_5: extra simulations with the active axon stubs (that were not included in the thesis work). Only symmetrical waveforms were tested.
- within the subfolders the experiments 2, 3 and 4 are structured as follows:
  - bkg: 100 ms simulation without extracellular stimulation
  - init_state: 1.5 s simulation without extracellular stimulation; the extracellular stimulation tests will always start from this initial state
  - pattern_0: single layer stimulation layer 4, symmetrical waveform
  - pattern_4: single layer stimulation layer 4, asymmetrical waveform
  - pattern_5: multi-layer stimulation 1 return electrorde layer 4 to 2/3, symmetrical waveform
  - pattern_6: multilayer stimulation 1 return electrode layer 4 to 2/3, asymmetrical waveform
  - pattern_7: multilayer stimulation 2 return electrodes layer 4 to 2/3, symmetrical waveform
  - pattern_8: multilayer stimulation 2 return electrodes layer 4 to 2/3, asymmetrical waveform
- the subfolders of experiment 5 are explained in a separate readme within the folder. All applied waveforms are symmetrical waveforms. Figure electrode_numbers_and_patt_explained can also help to understand the numbering of electrodes and patterns in exp 5.

### 4. Toolbox
- contains all extra helperfunctions; here follows an overview
- analyse functions and data_analyser: code that constructs the different plots:
- build_network_no_v1_edges: the build file used to build the networks for exp3: no virtual connections between the different neurons
- calcium_experiments: this code is a copy of the code applied in the calcium experiments and served as inspiration for the analysis code
- file_helper: contains helper functions to construct config files and retrieve parameters automatically
- hdf5 and hdf5_edge: contains functions to extract information from .h5 files
- hyperstim-db: code from the hyperstim project that served as inspiration for the analysis code.
- mask: code that loads the mask and applies it when building the network. This mask can be found in components/stimulation/comsol and is a simulation in comsol were all electrodes have a potential of 1; this way the exact position of the electrodes can be defined. With the mask.py code is it made sure no neurons coincide with the location of the electrodes.
- spikes_helper: helper functions that retrieves spikerates from the sonata files
- waveform.py: code to construct new waveforms

### 5. Virtual mice mask
- contains 3 instances (each generated with different random seed) of SONATA network files describing the V1 model. 
- these 'mice'/'networks' were used in experiment 2 and 4

### 6. Virtual mice mask no v1 edges
- contains 3 instances (each generated with different random seed) of SONATA network files describing the V1 model, without virtual connections. 
- these 'mice'/'networks' were used in experiment 3

### 7. Build_network.py
- use to rebuild networks
- in contrast to the original build file in bmtk, this file adjusts the node positions of neurons making sure they cannot coincide with the location of the probes

### 8. Run_bionet.py files
- run_bionet.py - use to run simulations
- run_bionet_axon.py - use to run simulations with active axon stubs. Properties of the axons as well as number of stubs can be adjusted in this file

## Insert active axon properties
To insert extra axon stubs and give them extra properties modify the run_bionet_axon.py file. The two most important functions applied are:

- The function fix_axon_peri_multiple_stubs allows to specify the number of stubs and their diameter

- The function set_params_peri_axon allows for inserting the equations of the new modfile (that needs to be added to components/mechanisms/modfiles and compiled) and changing the conductances, axial resistance values and membrane capacitance.

## Running a simulation

The cell models use some customized NEURON channels that needs to be compiled on the HPC using the appropiate version of NEURON. This only needs to be done once - 
unless you plan to upgrade to a different version of NEURON. 
```
cd biophys_components/mechanisms
nrnivmodl modfiles 
cd ../.. 
```
For most clusters the simulation can be started using the following command inside a SLURM/MOAB/etc. script (replace N with the number of cores to use)
``` 
mpirun -np N nrniv -mpi -python run_bionet.py config.json 
```

## Simulation output

By default once the simulation has started running it will create an "output/" folder containing the simulation results. (WARNING: if an existing "output/" folder from
a previous simulation already exists bmtk will overwrite it). 

"output/log.txt" will keep a running tally of the simulation and can be good to check-on (eg $ tail -f output/log.txt) while the simulation is running. The log text is also printed in the terminal. When completed, the spike trains for all the V1 cells will be stored in "output/spikes.h5". The spikes are stored according to the SONATA format 
(https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md#spike-file), and can be read using tools like pysonata, libsonata, or any hdf5 
API (eg h5py).

You can change where and how the output is stored in the config.json files under the "outputs" section.


## Modifying the simulation

The conditions under which the simulation is ran is set in config.json files.

Information about the configuration files can be found in the [SONATA documentation](https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md#tying-it-all-together---the-networkcircuit-config-file).
Also see [here](https://github.com/AllenInstitute/sonata/tree/master/examples) and [here](https://github.com/AllenInstitute/bmtk/tree/develop/docs/examples) for various
examples of different types of simulations.


## V1 Stimuli

### Original inputs

The V1 model is being stimulated by feed-forward synaptic inputs from two population of cells called LGN and BKG. The LGN cells represent a model of the LGN in the thalamus, while BKG cells represent generic background inputs from higher cortical areas. The V1 simulations use existing spike-train files for both sets of cells, which are then used to determine the timing synaptic activity for the post-synaptic V1 cells. 

The LGN and BKG spike-train files are saved in ../lgn_stimulus/results/* and ../bkg_inputs/results/* (in the [dropbox]((https://www.dropbox.com/sh/w5u31m3hq6u2x5m/AACpYpeWnm6s_qJDpmgrYgP7a?dl=0))), these folders also contains scripts and instructions for how to generate new stimulus for both sets of cells. Modify which spike-train files are read can be done by editing the "inputs" sections of config.json

### Extracellular potentials

In addition to the natural inputs coming from the LGN and the rest of the brain, an artificial input in the form of extracellular stimulation can be implemented with the comsol module in BioNet. To do this, you need a .txt file containing the FEM solution as well as a .csv file containing the waveform of the pulses. For more information about this, refer to the [comsol module](../bmtk/simulator/bionet/modules/comsol.py) and [bio_comsol example](../examples/comsol/README.md).


## Rebuilding the V1 Models

The Biophysically detailed recurrent V1 is already saved in the "Biophysical_network" folder, in the latest SONATA format. Should one want to rebuild the network 
files, the folder "network_builder" contains scripts and metadata to do so. The script to do so is called build_network.py, and requires Python 3.6+ and BMTK 
installed. On a single core machine the process of rebuilding the entire network from scratch can take a few hours at minimum. However one can also rebuild only 
parts of the network (ex regenerating the connection matrix but leaving the cell locations the same), or even build a smaller subsampled version of the network. 
See the README.txt inside network_builder folder for more instructions.

To rebuild the network:
```
$ python build_network.py
```

This can also be done in parallel if using an HPC with MPI installed
```
$ mpirun -np N python build_network.py
```

WARNING: RUNNING build_network.py WILL OVERWRITE EXISTING NETWORK FILES IN THE network/ FOLDER.


The following describes some of the base rules for how the network was built:

### Setting up a high-level description of populations
	
The populations in [build_files/biophys_props/v1_node_models.json](build_files/biophys_props/v1_node_models.json) are defined as e.g.:
```
    "i4Sst":   {"ncells":2384,
                "labels": {"ei":"i", "location":"VisL4"},
                "depth_range":[310,430],
                "criteria":{"ei":"i", "location":"VisL4",   
                "cre_line":["Sst"], "express":True}},
```
where "criteria" are used to choose models from the available models in the bio_models_prop.csv. The number of cells in each population is proportional to our 
expectation about the relative occurrence of each cre-line in mouse V1. The numbers of cells are pre-calculated in the file column_structure.xls. For the domain we 
use a right circular cylinder with a total radius of 845 um. The core cylinder within 400 um is populated with the biophysical cell models and the remaining annulus 
with the lif models.

### Assign somatic coordinates

Cells for each population are uniformly distributed within a cylindrical domain (400 µm radius core and 845 µm radius annulus) and within the specified "depth_range". Via the mask function, we avoid neurons to be placed at the location of the electrodes.


### Assign rotation_angle_yaxis

Uniformly draw rotational_angle_yaxis for each cell.


### Assign cell models
	
Assign model to cells based on their y_soma position. A model may be assigned to a particular cell if that model's morphology does not significantly stick out of 
the pia when placed at the cell's somatic location. Currently, we allow for dendrites to stick out of the pia with a tolerance = 100 um because there are not enough 
short L23 excitatory cells. Determine the acceptable yrange = [ymin,ymax] for each model based on its dendritic span [delta_ymin,delta_ymax]. 
	
For a given cell position, find models for which y_soma is in yrange. Randomly choose from the allowed models with a Gaussian 
(prob = np.exp(-(y_soma-y_spec)**2/length_sigma**2/2)) probability density function (length_sigma = 20 um). The probability density function is chosen such that we 
preferentially choose a model, which has a y_spec closest to the y_soma of a particular cell. For the periphery we substitute the biophysical models with the lif models 
according to the mapping bio2lif_mapping.csv


