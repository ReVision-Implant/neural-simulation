#!/bin/bash

# Calling this script in the terminal will build all the networks specified in the lines starting with mpirun
# $ bash /path/to/build.sh

module load mpi
cd "${0%/*}"        # Change the working dir to where this script is located

mpirun -np 12 python build_network.py --fraction 0.001 -o networks_25/network0 --rng-seed 100  # !!! Change the arguments according to your needs !!!
# You can add as many lines as you want, e.g.
# mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network1 --rng-seed 101
# mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network2 --rng-seed 102
# mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network3 --rng-seed 103
# mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network4 --rng-seed 104


