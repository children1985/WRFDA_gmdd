#!/bin/csh

### Project name
#PBS -A NMMM0015

### Job name
#PBS -N tl_test

### Wallclock time
#PBS -l walltime=00:10:00

### Queue
#PBS -q share

### Merge output and error files
#PBS -j oe                    

### Select 2 nodes with 36 CPUs, for 72 MPI processes 
#PBS -l select=1:mpiprocs=1


./linear_check.exe >& log.log

exit
