#!/bin/bash
source /users/students/r0754386/Documents/bmtk/bin/activate
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1
mpirun -np 12 nrniv -mpi -python run_bionet.py exp3-/config/1/-/10/config_01-_10_1.json



