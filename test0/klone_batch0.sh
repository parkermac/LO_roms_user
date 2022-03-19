#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=1
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=4

## Walltime 
#SBATCH --time=05:00:00

## Memory per node
#SBATCH --mem=100G

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort/
export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

mpirun -np 4 /mmfs1/gscratch/macc/parker/LO_roms_user/test0/romsM /mmfs1/gscratch/macc/parker/LO_roms_user/test0/roms_upwelling.in > /mmfs1/gscratch/macc/parker/LO_roms_user/test0/roms_log.txt

