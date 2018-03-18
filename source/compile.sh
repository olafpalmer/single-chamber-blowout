#!/bin/bash


gfortran -c -static io_mod.f90
gfortran -c -static single_chamber.f90
gfortran -static blowout_calculation.f90 io_mod.o single_chamber.o -o executable 
