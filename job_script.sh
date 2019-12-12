#!/bin/bash

# number of tasks = number of Processes = 1 for openmp
#SBATCH --ntasks 1

# amount of memory per core in MB
#SBATCH --mem-per-cpu=1200

# Real Time the job should run (HH:MM:SS)
#SBATCH --time 00:30:00

# number of cores per tasks = number of (openmp) threads per process
#SBATCH --cpus-per-task=4

# Name of the Job in the queue (arbitrary choice)
#SBATCH --job-name Hello_World

# where the output of the job goes (filename)
#SBATCH --output job.out

# where stderr goes (filename)
#SBATCH --error job.err

# Use the reservation of our course:
#SBATCH -A kurs00035
#SBATCH -p kurs00035
#SBATCH --reservation=kurs00035

# see man sbatch for more parameters of jobscripts

# here the jobscript starts
# the following commands will be executed on the compute nodes

#reset environment
module purge
# load the compiler needed
module load gcc

# set the number of openmp hreads according to the specification above
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# your executable here
srun ./a.out
