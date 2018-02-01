#!/bin/sh
basn=`basename $1 .cu`
nvcc -I/usr/local/cuda-9.1/include $1 -L /usr/local/cuda-9.1/lib64/ -lcublas -lcusparse -o $basn.x
