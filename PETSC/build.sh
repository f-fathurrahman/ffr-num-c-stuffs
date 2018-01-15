#!/bin/bash
basn=`basename $1 .c`

# Using static lib
mpicc -Wall $1 -I/usr/local/petsc-3.8.3_openmpi_gnu/include \
/usr/local/petsc-3.8.3_openmpi_gnu/lib/libpetsc.a \
-lblas -llapack -lm -lX11 -ldl -o $basn.x

