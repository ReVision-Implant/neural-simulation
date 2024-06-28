# The single neuron model

The working of the single neuron model is very similar to that of the V1 model. Therefore, the README file here is quite short. For details refer to the v1_Anke on the main branch. 

The comsol files that calculated the extracellular potentials can be found under components/stimulation/comsol/electrode_x_small. !caveat: the comsol set-up here is different from the v1 model (much simpler, two rows of 4 electrodes, no probes taken into account). Mostly only electrode 0 is used here.

For the single neuron stimulation we work in 3 steps:
- First the network is build with the builder.py. Here specify the location of your neurons wrt to the electrodes. The position of the electrodes in the single neuron model can be found under components/stimulation/electrodes. Also the specific celltype (active properties of the soma) can be chosen and changed here. The analyser file contains code to plot the location of the electrode versus the 2 neurons. Mostly we only look at one neuron and ignore the other one.
- Secondly the simulation is set up with the simulator.py file. Here specify the folder/name of the simulation, the network and duration of the simulation. Furthermore, here we indicate that we want to track the membrane potential as well. If you need your own modfiles for specifying the equations of the ionchannels you will need to copy them to the newly created components folder and compile them.
- Third: run the simulation with the run.py file. Here specify the number of stubs for the axon and the properties that need to be added.

## Active axon properties
To insert extra axon stubs and give them extra properties modify the run.py file. The two most important functions applied are:

    The function fix_axon_peri_multiple_stubs allows to specify the number of stubs and their diameter

    The function set_params_peri_simpl_hh allows for inserting the equations of the new modfile (that needs to be added to components/mechanisms/modfiles and compiled) and changing the conductances, axial resistance values and membrane capacitance.
