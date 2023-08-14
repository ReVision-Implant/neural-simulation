# Experiment 0

This folder contains the files corresponding to the experiment with monopolar stimulation.

## File structure

- config/ - contains the configuration files. Each defines the parameters of a single simulation. Structured like electrode>stimulation_type>amplitude>config.json
- output/ - contains the output (config, log, spikes) of the simulations. Structured identically as config/.
- config_creator.py - Script used to generate config files so you don't have to manually change parameters in tens of files.
- plot_output.py - Used for plotting the results.
- runall.sh - Bash script used to run all simulations (that haven't been run yet). Running this in the command line allows you to run multiple simulations sequentially.
- testall.sh - Bash script that prints out how many simulations have already been done.
