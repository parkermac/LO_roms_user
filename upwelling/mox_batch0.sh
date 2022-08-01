#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=4

## Walltime 
#SBATCH --time=00:20:00

## Memory per node
#SBATCH --mem=128G

module load icc_17-impi_2017
module load netcdf_fortran+c_4.4.1.1-icc_17

RUN_DIR=/gscratch/macc/parker/LO_roms_user/upwelling
mpirun -np 4 $RUN_DIR/romsM $RUN_DIR/roms_upwelling.in > $RUN_DIR/roms_log.txt

