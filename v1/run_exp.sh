#!/bin/bash

# Calling this script in the terminal will run or test all the config files in the specified experiment
# $ bash run_exp.sh [arg1] [arg2]
#   arg1: name of experiment folder, e.g. 'exp0'
#   arg2: 'run' or 'test'
#       'run' will run simulations that have not been ran before
#       'test' will print the number of completed and total simulations

# Initialisation
module load mpi
cd "${0%/*}"        # Change the working dir to where this script is located       
echo '_______________'
TOTAL=0
DONE=0

# Find all files in [arg1]/config/ (and subdirectories) and loop over the file names, checking each time if output already exists.
config_path="$(pwd)/$1/config/"                             # Navigate to v1/[arg1]
for config_file in $(find $config_path -type f -print)      # Loop over all files in config_path and its subdirectories
do
    let TOTAL++                                             # Increase variable by one
    output_path=${config_file/config/"output"}              # Replace "config" by "output" in the path
    file_name="$(basename "$output_path")"                  # Get the final component of the output path (everything after the last /)
    file_name="${file_name:7:-5}"                           # Remove "config" and ".json" from file_name
    output_path="$(dirname "$output_path")"                 # Get parent dir of output_path
    output_path="${output_path}/${file_name}/spikes.csv"    # Concatenate
    if [[ -f $output_path ]]; then                          # If there are already files in the output path, then
        if [ $2 == "run" ]; then                            # If second input argument is 'run', then
            echo "$file_name already exists"                # Print file name
        fi
        let DONE++                                          # Increase variable by one
    else 
        echo "$file_name does not exist"                    # Print file name
        if [ $2 == "run" ]; then                            # If second input argument is 'run', then
            nice -10 mpirun -np 1 nrniv -mpi -python run_bionet.py $config_file       # Run simulation
        fi
    fi
done
if [ $2 == "test" ]; then                                   # If second input argument is 'test', then
    echo "$DONE/$TOTAL"                                     # Print number of completed and total simulations.
fi