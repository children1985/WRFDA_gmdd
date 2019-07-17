#!/bin/csh

set echo

set nclib = $NETCDF/lib
set ncinc = $NETCDF/include

set fc = ifort
set fflg1 = "-convert little_endian -qopenmp"
set fflg2 = "-I$ncinc -L$nclib -lnetcdf -lm"

set fflg = ( $fflg1 $fflg2 )

ifort $fflg -c wrfout.f90

ifort -o wrfout.exe $fflg wrfout.o

