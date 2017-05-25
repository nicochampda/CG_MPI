#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 8:00
#SBATCH -n 16

module load gcc openmpi
mpirun ./conjgrad 256
mpirun ./conjgrad 384
mpirun ./conjgrad 512
mpirun ./conjgrad 640
mpirun ./conjgrad 768
mpirun ./conjgrad 896
mpirun ./conjgrad 1024
mpirun ./conjgrad 1152
mpirun ./conjgrad 1280
mpirun ./conjgrad 1408
