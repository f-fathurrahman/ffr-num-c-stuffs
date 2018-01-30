export SLEPC_DIR=/usr/local/slepc-3.8.2_openmpi_gnu/
export PETSC_DIR=/usr/local/petsc-3.8.3_openmpi_gnu/

export INC_SLEPC=$SLEPC_DIR/include/
export INC_PETSC=$PETSC_DIR/include/

export LIB_SLEPC="-L$SLEPC_DIR/lib -lslepc"
export LIB_PETSC="-L$PETSC_DIR -lpetsc"

basnam=`basename $1 .c`
mpicc -I$INC_SLEPC -I$INC_PETSC $1 $LIB_SLEPC $LIB_PETSC -o $basnam.x

