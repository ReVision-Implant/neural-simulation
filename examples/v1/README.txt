Biophysical_network
===================

Contains the files necessary to run a full simulation of the biophysically detailed V1 network with thalamacortical (LGN) and background (BKG) inputs. Is designed
to run using the Allen Institute [Brain Modeling Toolkit](https://github.com/AllenInstitute/bmtk); but the network, input and config files are in the 
[SONATA data format](https://github.com/AllenInstitute/sonata) for use with other tools that support SONATA.


Requirements
------------
Python 3.6+
NEURON 7.4+
BMTK 0.0.8+

At the moment it also requires access to a HPC cluster, as the network is too large to instantiate in even the most powerful of desktop machines (by 2019 standards). 
A cluster with 100 Xeon cores will typically take a few hours to run through the entire 3.0 second simulation. 


File Structure
--------------
* network/ - Contains SONATA network files describing the full V1 model
* network_precomputed_syns/ - Similar to above but with the synaptic locations and synaptic weights fully specified. In the network/*edges_types.csv the columns weight_function,
     target_sections and distance_range are used by bmtk to calculate where and how to place synapses between every source/target pair. In these files the synaptic locations
     and weights are stored in the edges.h5 files, more complient with the SONATA format.
* components/ - parameters, morphology, model and mechanisms files used to instantiate individual cells and synapses
* build_files/ - scripts and properties for rebuilding the V1 network.


Running the simulation
----------------------
The cell models use some customized NEURON channels that needs to be compiled on the HPC using the appropiate version of NEURON. This only needs to be done once - 
unless you plan to upgrade to a different version of NEURON.
  $ cd biophys_components/mechanisms
  $ nrnivmodl modfiles
  $ cd ../..

For most clusters the simulation can be started using the following command inside a SLURM/MOAB/etc. script (replace N with the number of cores to use)
  $ mpirun -np N nrniv -mpi -python run_bionet.py config.json


Simulation output
-----------------
By default once the simulation has started running it will create an "output/" folder containing the simulation results. (WARNING: if an existing "output/" folder from
a previous simulation already exists bmtk will overwrite it). 

"output/log.txt" will keep a running tally of the simulation and can be good to check-on (eg $ tail -f output/log.txt) while the simulation is running. When completed
the spike trains for all the V1 cells will be stored in "output/spikes.h5". The spikes are stored according to the SONATA format 
(https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md#spike-file), and can be read using tools like pysonata, libsonata, or any hdf5 
API (eg h5py).

You can change where and how the output is stored in "simulation_config.json" under the "outputs" section.


Modifying the simulation
------------------------
The conditions under which the simulation is ran is set in two json configuration files:
* circuit_config.json - contains information required to instaniate the network; what set(s) of cells and synapses to use, location to cell and synaptic parameters, etc.
* simulation_config.json - information about the simulation; run time, dt, inputs and what gets recorded.

Information about the configuration files can be found in the [SONATA documentation](https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md#tying-it-all-together---the-networkcircuit-config-file).
Also see [here](https://github.com/AllenInstitute/sonata/tree/master/examples) and [here](https://github.com/AllenInstitute/bmtk/tree/develop/docs/examples) for various
examples of different types of simulations.


V1 Stimuli
----------
The V1 model is being stimulated by feed-forward synaptic inputs from two population of cells called LGN and BKG. The LGN cells represent a model of the LGN in the 
thalamus, while BKG cells represent generic background inputs from higher cortical areas. The V1 simulations use existing spike-train files for both sets of cells,
which are then used to determine the timing synaptic activity for the post-synaptic V1 cells.

The LGN and BKG spike-train files are saved in ../lgn_stimulus/results/* and ../bkg_inputs/results/*, these folders also contains scripts and instructions for how
to generate new stimulus for both sets of cells. Modify which spike-train files are read can be done by editing the "inputs" sections of config.json



Rebuilding the V1 Models
---------------
The Biophysically detailed recurrent V1 is already saved in the "Biophysical_network" folder, in the latest SONATA format. Should one want to rebuild the network 
files, the folder "network_builder" contains scripts and metadata to do so. The script to do so is called build_network.py, and requires Python 3.6+ and BMTK 
installed. On a single core machine the process of rebuilding the entire network from scratch can take a few hours at minimum. However one can also rebuild only 
parts of the network (ex regenerating the connection matrix but leaving the cell locations the same), or even build a smaller subsampled version of the network. 
See the README.txt inside network_builder folder for more instructions.

To rebuild the network:
  $ python build_network.py

This can also be done in parallel if using an HPC with MPI installed
  $ mpirun -np 2 python build_network.py

WARNING: RUNNING build_network.py WILL OVERWRITE EXISTING NETWORK FILES IN THE network/ FOLDER.


The following describes some of the base rules for how the network was built:

### Setting up a high-level description of populationsx xls
	The populations in build_files/biophys_props/pop_queries.json are defined as e.g.:

        "i4Sst":   {"ncells":2384,
                    "labels": {"ei":"i", "location":"VisL4"},
                    "depth_range":[310,430],
                    "criteria":{"ei":"i", "location":"VisL4", "cre_line":["Sst"], "express":True}},

where "criteria" are used to choose models from the available models in the bio_models_prop.csv. The number of cells in each population is proportional to our 
expectation about the relative occurrence of each cre-line in mouse V1. The numbers of cells are pre-calculated in the file column_structure.xls. For the domain we 
use a right circular cylinder with a total radius of 845 um. The core cylinder within 400 um is populated with the biophysical cell models and the remaining annulus 
with the lif models.

### Assign somatic coordinates

Cells for each population are uniformly distributed within a cylindrical domain and within the specified "depth_range".


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



### Test and visualize the node construction

* Build a network with one cell per model, construct and save segment coordinates of each morphology to build_model/cells_peri/biophysical/morph_segs
* analyze/plots the depth distributions of the models and their morphologies with build_model/V1/plot_model_distr.ipynb


 ### Information

 networks_rebuilt contains network only including 0.10 fraction of neurons.
 The y-component represents the axial direction of the column, representing L1->6 with increasing y


ssh -J r0754386@ssh.esat.kuleuven.be r0754386@eridani2.esat.kuleuven.be