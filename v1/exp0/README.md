# Experiment 0

This folder contains the files corresponding to the experiment with monopolar stimulation.

## File structure

- config/ - contains the `configuration.json` files. Structured like `electrode>stimulation_type>amplitude>config.json`. Each `.json` file defines the parameters of a single simulation. 
- output/ - Structured identically as config/. contains the output folder containing 
    - copy of the configuration file
    - log.txt
    - spikes 
    - any other output variables you might have specified
- config_creator.py - Used to generate config files so you don't have to manually change parameters in tens or hundreds of files.
- plot_output.py - Used for plotting the results.
