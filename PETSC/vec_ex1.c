#include <petscvec.h>

int main( int argc, char** argv )
{
  PetscErrorCode ierr;
  ierr = PetscInitialize( &argc, &argv, NULL, NULL ); CHKERRQ( ierr );

  PetscPrintf( PETSC_COMM_WORLD, "Program started\n" );

  int Nprocs;
  ierr = MPI_Comm_size( PETSC_COMM_WORLD, &Nprocs );
  PetscPrintf( PETSC_COMM_WORLD, "Program is using %d procs\n", Nprocs );

  // Create a vector (parallel)
  Vec x;
  PetscInt N = 10;
  VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, N, &x );

  // Duplicate vector
  Vec y, w;
  VecDuplicate( x, &y );
  VecDuplicate( x, &w );

  // array of vectors, by duplicating
  Vec *z;
  VecDuplicateVecs( x, 3, &z );

  // set their values
  PetscScalar ONE = 1.0, TWO = 2.0, THREE = 3.0;
  VecSet( x, ONE );
  VecSet( y, TWO );
  VecSet( z[0], ONE );
  VecSetRandom( z[1], NULL );
  VecSet( z[2], 2.5 );

  if( N <= 10 ) {
    VecView( x, PETSC_VIEWER_STDOUT_WORLD );
    VecView( y, PETSC_VIEWER_STDOUT_WORLD );
    VecView( z[0], PETSC_VIEWER_STDOUT_WORLD );
    VecView( z[1], PETSC_VIEWER_STDOUT_WORLD );
    VecView( z[2], PETSC_VIEWER_STDOUT_WORLD );
  }

  // Free memory
  VecDestroy( &x );
  VecDestroy( &y );
  VecDestroy( &w );
  VecDestroyVecs( 3, &z );


  PetscPrintf( PETSC_COMM_WORLD, "Program will end soon\n");

  ierr = PetscFinalize(); CHKERRQ( ierr );
  return 0;
}
