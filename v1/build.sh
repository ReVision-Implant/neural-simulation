#!/bin/bash
source /users/students/r0754386/Documents/bmtk/bin/activate
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1
mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network0 --rng-seed 100
mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network1 --rng-seed 101
mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network2 --rng-seed 102
mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network3 --rng-seed 103
mpirun -np 12 python build_network.py --fraction 0.25 -o networks_25/network4 --rng-seed 104


