#!/bin/bash

# Initialisation
source /users/students/r0754386/Documents/bmtk/bin/activate # Activate environment
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1 # So the script can be run from any dir
echo '_______________'

# Find all files in the config/ dir and loop over the file names, checking each time if output already exists.
for config_file in $(find /users/students/r0754386/Documents/bmtk/examples/v1/exp1/config -type f -print)
do
    output_path=${config_file/config/"output"}      # Replace "config" to "output" in the path
    file_name="$(basename "$output_path")"                  # Get the final component of the output path (everything after the last /)
    file_name="${file_name:7:-5}"                           # Remove "config" and ".json" from file_name
    output_path="$(dirname "$output_path")"                 # Get parent dir of output_path
    output_path="${output_path}/${file_name}/spikes.csv"    # Concatenate
    if [[ -f $output_path ]]; then                          # If there are already files in the output path
        echo "$output_path already exists"
    else 
        echo "$output_path does not exist"
        nice -10 mpirun -np 12 nrniv -mpi -python run_bionet.py $config_file       # Run simulation
    fi
done