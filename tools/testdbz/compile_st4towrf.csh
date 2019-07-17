#!/bin/csh

set echo

set nclib = $NETCDF/lib
set ncinc = $NETCDF/include

set fc = ifort
set fflg1 = "-convert little_endian"
set fflg2 = "-I$ncinc -L$nclib -lnetcdf -lm"

set fflg = ( $fflg1 $fflg2 )

ifort $fflg -c st4towrf.f90

ifort -o st4towrf.exe $fflg st4towrf.o

