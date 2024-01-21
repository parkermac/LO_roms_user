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

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort
NCDIR=/gscratch/macc/local/netcdf-icc
export LD_LIBRARY_PATH=${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}

RUN_DIR=/gscratch/macc/parker/LO_roms_user/upwelling
mpirun -np 4 $RUN_DIR/romsM $RUN_DIR/roms_upwelling.in > $RUN_DIR/roms_log.txt

