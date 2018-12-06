#!/bin/bash

# Modify as you need
export filepath=/net/scratch3/guofan/trinity/turbulent-sheet3D-mixing-sigma100/tracer
export particle=electron
tstep=0

# May not need to change
mpi_size=1
fpath=$filepath/T.$tstep
input_file=$fpath/${particle}_tracer_sorted.h5p

cd ../

srun -n $mpi_size ./check_hdf5 $input_file /Step#$tstep

cd config
