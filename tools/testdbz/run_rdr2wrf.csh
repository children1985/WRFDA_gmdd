#!/bin/csh

### Project name
#PBS -A NMMM0015

### Job name
#PBS -N rdr2wrf

### Wallclock time
#PBS -l walltime=00:30:00

### Queue
#PBS -q share

### Merge output and error files
#PBS -j oe                    

### Select 2 nodes with 36 CPUs, for 72 MPI processes 
#PBS -l select=1:mpiprocs=1

./rdr2wrf.exe >& rdr2wrf.log

exit
