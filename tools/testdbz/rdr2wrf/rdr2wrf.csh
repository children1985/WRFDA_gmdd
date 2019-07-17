#!/bin/csh

### Project name
#PBS -A NMMM0015

### Job name
#PBS -N rdr2wrf

### Wallclock time
#PBS -l walltime=00:10:00

### Queue
#PBS -q regular

### Merge output and error files
#PBS -j oe                    

### Select 2 nodes with 36 CPUs, for 72 MPI processes 
#PBS -l select=1:ncpus=36:ompthreads=36

set echo

#set nodes = 1
#set cores = 6
#@ np = $nodes * $cores

#cp namelist.input.temp namelist.input

cp config.dat.ob config.dat
./cal_wrf_dbz.exe

cp config.dat.fg config.dat
./cal_wrf_dbz.exe

cp config.dat.an config.dat
./cal_wrf_dbz.exe

exit
