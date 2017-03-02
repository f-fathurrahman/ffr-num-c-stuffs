export INC_SLEPC=/home/efefer/mysoftwares/slepc-3.7.3/include
export INC_PETSC=/home/efefer/mysoftwares/petsc-3.7.5/include
export LIB_SLEPC="-L/home/efefer/mysoftwares/slepc-3.7.3/lib -lslepc"
export LIB_PETSC="-L /home/efefer/mysoftwares/petsc-3.7.5/lib -lpetsc"

basnam=`basename $1 .c`
mpicc -I$INC_SLEPC -I$INC_PETSC $1 $LIB_SLEPC $LIB_PETSC -o $basnam.x

