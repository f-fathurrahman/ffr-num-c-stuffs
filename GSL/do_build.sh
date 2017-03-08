export LIB_GSL="-lgsl -lgslcblas -lm"

basnam=`basename $1 .c`
gcc $1 $LIB_GSL -o $basnam.x

