#!/bin/bash
basn=`basename $1 .c`

mpicc -Wall -I/usr/local/petsc-3.8.3_openmpi_gnu/include $1 \
/usr/local/petsc-3.8.3_openmpi_gnu/include/libpetsc.a -o $basn.x

