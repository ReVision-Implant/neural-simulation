#!/bin/bash

# Calling this script in the terminal will run all config files specified in the lines starting with mpirun
# $ bash /path/to/run.sh

module load mpi
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # Change the working dir to where this script is located         
mpirun -np 12 nrniv -mpi -python run_bionet.py exp0/config/config0.json         # !!! Change the last argument to the path of the config file you want to simulate !!!
# You can add as many lines as you want, e.g.
# mpirun -np 12 nrniv -mpi -python run_bionet.py exp0/config/config1.json
# mpirun -np 12 nrniv -mpi -python run_bionet.py exp0/config/config2.json
# mpirun -np 12 nrniv -mpi -python run_bionet.py exp0/config/config3.json


