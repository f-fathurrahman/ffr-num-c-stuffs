#include <petsc.h>

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


  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  //
  // Dot products: VecDot and VecMDot
  //
  PetscScalar dot;
  VecDot( x, y, &dot );
  PetscPrintf( PETSC_COMM_WORLD, "dot = %18.10f\n", dot );
  //
  PetscScalar dots[3];
  VecMDot( x, 3, z, dots );
  int i;
  for( i = 0; i < 3; i++ ) {
    PetscPrintf( PETSC_COMM_WORLD, "dots[%d] = %18.10f\n", i, dots[i] );
  }

  
  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  PetscPrintf( PETSC_COMM_WORLD, "Vector length: %d\n", N );

  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  // 
  // Maximum and minimum value
  //
  PetscReal maxval, minval;
  PetscInt maxidx, minidx;
  //
  VecMax( z[1], &maxidx, &maxval );
  PetscPrintf( PETSC_COMM_WORLD, "Vec z[1]: maxidx = %d, maxval = %18.10f\n", maxidx, maxval );
  //
  VecMin( z[1], &minidx, &minval );
  PetscPrintf( PETSC_COMM_WORLD, "Vec z[1]: minidx = %d, minval = %18.10f\n", minidx, minval );

  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  //
  // Scale and norm
  //
  VecScale( x, 1.1 );
  PetscReal norm;
  VecNorm( x, NORM_2, &norm );
  PetscPrintf( PETSC_COMM_WORLD, "Vec x: after scale by 1.1 norm = %18.10f\n", norm );
  PetscReal norm0;
  norm0 = PetscSqrtReal( 1.1*1.1*N );
  PetscPrintf( PETSC_COMM_WORLD, "Norm should be %18.10f\n", norm0 );

  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  PetscPrintf( PETSC_COMM_WORLD, "PETSC_SMALL = %18.10e\n", PETSC_SMALL );


  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  //
  // Copy
  //
  VecCopy( x, w );
  PetscReal normw;
  VecNorm( w, NORM_2, &normw );
  PetscPrintf( PETSC_COMM_WORLD, "Vec w: copied from x, normw = %18.10f\n", normw );
  PetscReal normw0;
  normw0 = PetscSqrtReal( 1.1*1.1*N );
  PetscPrintf( PETSC_COMM_WORLD, "normw should be %18.10f\n", normw0 );


  PetscReal v;

  VecAXPY( y, THREE, x );
  VecNorm( y, NORM_2, &norm );
  PetscPrintf( PETSC_COMM_WORLD, "Vec w: copied from x, normw = %18.10f\n", normw );
  v    = norm - 8.0*PetscSqrtReal((PetscReal)N); if (v > -PETSC_SMALL && v < PETSC_SMALL) v = 0.0;
  PetscPrintf(PETSC_COMM_WORLD,"VecAXPY %g\n", (double)v);

  VecAYPX( y,TWO, x );
  VecNorm( y, NORM_2, &norm);
  v    = norm-18.0*PetscSqrtReal((PetscReal)N);
  PetscPrintf(PETSC_COMM_WORLD,"VecAYPX %g\n", (double)v );

  VecSwap( x, y );
  VecNorm( y, NORM_2, &norm );
  v    = norm-2.0*PetscSqrtReal( (PetscReal)N );
  PetscPrintf( PETSC_COMM_WORLD, "VecSwap  %g\n", (double)v );
  VecNorm( x, NORM_2, &norm );
  v = norm - 18.0*PetscSqrtReal( (PetscReal)N );
  PetscPrintf( PETSC_COMM_WORLD, "VecSwap  %g\n", (double)v );

  VecWAXPY( w, TWO, x, y );
  VecNorm( w, NORM_2, &norm );
  v    = norm - 38.0*PetscSqrtReal( (PetscReal)N );
  PetscPrintf( PETSC_COMM_WORLD, "VecWAXPY %g\n", (double)v);

  VecPointwiseMult( w, y, x );
  VecNorm( w, NORM_2, &norm );
  v    = norm - 36.0*PetscSqrtReal( (PetscReal)N );
  PetscPrintf(PETSC_COMM_WORLD,"VecPointwiseMult %g\n",(double)v);

  VecPointwiseDivide(w,x,y);
  VecNorm(w,NORM_2,&norm);
  v    = norm-9.0*PetscSqrtReal((PetscReal)N);
  PetscPrintf(PETSC_COMM_WORLD,"VecPointwiseDivide %g\n",(double)v);

  //
  // Free memory
  //
  VecDestroy( &x );
  VecDestroy( &y );
  VecDestroy( &w );
  VecDestroyVecs( 3, &z );


  PetscPrintf( PETSC_COMM_WORLD, "\n" );
  PetscPrintf( PETSC_COMM_WORLD, "Program will end soon\n");

  ierr = PetscFinalize(); CHKERRQ( ierr );
  return 0;
}
