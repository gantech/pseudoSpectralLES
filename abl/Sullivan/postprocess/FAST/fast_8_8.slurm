#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run multiple mpi
# applications within one batch job on TACC's
# Stampede system.
#----------------------------------------------------
#SBATCH -J FAST # Job name
#SBATCH -o FAST.o%j # Name of stdout output file(%j expands to jobId)
#SBATCH -e FAST.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p serial # Submit to the 'normal' or 'development' queue
#SBATCH -N 1 # Total number of nodes requested (16 cores/node)
#SBATCH -n 1 # Total number of mpi tasks requested
#SBATCH -t 12:00:00 # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=ganeshv@psu.edu
#SBATCH --mail-type=begin # email me when the job starts
#SBATCH --mail-type=end # email me when the job finishes# The next line is required if the user has more than one project
# #SBATCH -A A-yourproject # Allocation name to charge job against

./FASTfilesWrite 8 8 &> log.FASTfilesWrite_8_8

