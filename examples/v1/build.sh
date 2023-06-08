#!/bin/bash
source /users/students/r0754386/Documents/bmtk/bin/activate
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1
mpirun -np 12 python build_network.py --fraction 0.3183 -o square/network0 --rng-seed 100
mpirun -np 12 python build_network.py --fraction 0.3183 -o square/network1 --rng-seed 101
mpirun -np 12 python build_network.py --fraction 0.3183 -o square/network2 --rng-seed 102


