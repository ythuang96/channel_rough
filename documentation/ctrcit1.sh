#!/bin/bash
##      opciones del PBS; ver "man qsub" para mas informacion

# Name of job
#PBS -N ctrcittest1

# From JJ
#PBS -e out.err
##  PBS -o out.out

# Specify the shell types
#PBS -S /bin/bash

# Number of nodes (4 nodes with 2 CPUs each implies 16 here - the total number 
# of nodes passed to mpirun will be nodes@ppn)
#PBS -l nodes=16

# From JJ
#PBS -l mem=2gb

#PBS -l walltime=20:00:00

# Export all environment variables to the job
#PBS -V

# Specify the queue type
#PBS -q default

# Mail to the user when job terminates or aborts 
##PBS -m ae

#change the working directory (default is home directory)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# Write our some information on the job
echo Running on host `hostname`
echo Time is `date`

### Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

# Tell me which nodes it is run on
echo " "
echo This job runs on the following processors:
echo `cat $PBS_NODEFILE`
echo " "

#
# Run the mpi job
#

/opt/openmpi/bin/mpirun -np $NPROCS ROUGHv7


