export PETSC_DIR=/usr/local/petsc-3.8.3_openmpi_gnu

export INC_PETSC=$PETSC_DIR/include

export LIB_PETSC="-L${PETSC_DIR}/lib -lpetsc -lm"

echo $LIB_PETSC
basnam=`basename $1 .c`
mpicc -I$INC_PETSC $1 $LIB_PETSC -o $basnam.x

