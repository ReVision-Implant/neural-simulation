#!/bin/bash
source /users/students/r0754386/Documents/bmtk/bin/activate
module load mpi
cd /users/students/r0754386/Documents/bmtk/examples/v1
mpirun -np 12 nrniv -mpi -python run_bionet.py 9_z+/config_comsol_z+_gnd_+-.json
mpirun -np 12 nrniv -mpi -python run_bionet.py 9_z+/config_comsol_z+_gnd_+0.json
mpirun -np 12 nrniv -mpi -python run_bionet.py 9_z+/config_comsol_z+_ins_+-.json
mpirun -np 12 nrniv -mpi -python run_bionet.py 9_z+/config_comsol_z+_ins_+0.json



