#!/bin/csh

set echo
set fopt = "-convert big_endian"
set fopt = ""
#ifort $fopt -c oudbz.f90

#ifort $fopt -c linear_check.f90
#ifort -o linear_check.exe $fopt linear_check.o oudbz.o

#ifort $fopt -c linear_sensitivity.f90
#ifort -o linear_sensitivity.exe $fopt linear_sensitivity.o oudbz.o

ifort $fopt -c oudbz_test.f90
ifort $fopt -c qz_check.f90
ifort -o qz_check.exe $fopt qz_check.o oudbz_test.o
