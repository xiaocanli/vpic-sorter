#!/bin/bash

export filepath=/net/scratch3/guofan/trinity/turbulent-sheet3D-mixing-sigma100/tracer/

export particle=electron
nparticles=10000
ratio_emax=1
cd ../

mpirun -np 36 ./h5trajectory -d ${filepath} -o data/${particle}s_3.h5p \
    -n $nparticles -p ${particle} -r $ratio_emax

cd config
