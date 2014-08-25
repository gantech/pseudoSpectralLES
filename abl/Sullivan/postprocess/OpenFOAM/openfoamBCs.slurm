#!/bin/bash
#PBS -l size=12
#PBS -S /bin/bash
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -A TG-ATM110007


module swap PrgEnv-pgi PrgEnv-intel
cd $PBS_O_WORKDIR
./openFOAMBCsWrite