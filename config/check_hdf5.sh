#!/bin/bash

# Modify as you need
export filepath=/net/scratch3/xiaocanli/reconnection/open3d-full/tracer_test/
export particle=electron
tstep=13

# May not need to change
mpi_size=1
fpath=$filepath/T.$tstep
input_file=$fpath/${particle}_tracer_reduced_sorted.h5p

cd ../

srun -n $mpi_size ./check_hdf5 $input_file /Step#$tstep

cd config
