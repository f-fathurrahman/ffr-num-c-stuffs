#include <petscksp.h>

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";


int main( int argc, char** argv )
{
  PetscErrorCode ierr;

  ierr = PetscInitialize( &argc, &argv, (char*)0, help ); CHKERRQ( ierr );

  // Problem size
  PetscInt m = 8, n = 7;
  ierr = PetscOptionsGetInt( NULL, NULL, "-m", &m, NULL );
  ierr = PetscOptionsGetInt( NULL, NULL, "-n", &n, NULL );

  //
  // Matrix of linear system
  //
  Mat A;
  ierr = MatCreate( PETSC_COMM_WORLD, &A );
  ierr = MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, m*n, m*n );
  ierr = MatSetFromOptions( A );
  ierr = MatMPIAIJSetPreallocation( A, 5, NULL, 5, NULL );
  ierr = MatSeqAIJSetPreallocation( A, 5, NULL );
  ierr = MatSeqSBAIJSetPreallocation( A, 1, 5, NULL );
  ierr = MatMPISBAIJSetPreallocation( A, 1, 5, NULL, 5, NULL );

  //
  // Determined which rows of the matrix are locally owned
  //
  PetscInt Istart, Iend;
  ierr = MatGetOwnershipRange( A, &Istart, &Iend );
  printf("Istart = %d, Iend = %d\n", Istart, Iend );

  ierr = PetscFinalize(); CHKERRQ( ierr );

  printf("Program ended normally\n");
  return 0;
}
