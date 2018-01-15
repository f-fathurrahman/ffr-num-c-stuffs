#include <petscsys.h>

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";


int main( int argc, char** argv )
{
  PetscErrorCode ierr;
  Vec x, b, u;

  ierr = PetscInitialize( &argc, &argv, (char*)0, help );

  // Problem size
  PetscInt m = 8, n = 7;
  ierr = PetscOptionsGetInt( NULL, NULL, "-m", &m, NULL );
  ierr = PetscOptionsGetInt( NULL, NULL, "-n", &n, NULL );

  // Matrix of linear system
  Mat A;
  ierr = MatCreate( PETSC_COMM_WORLD, &A );
  ierr = MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, m*n, m*n );
  ierr = MatSetFromOptions( A );

}
