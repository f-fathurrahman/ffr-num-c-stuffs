#include <petsc.h>

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

  PetscLogStage stage;
  ierr = PetscLogStageRegister( "Assembly", &stage ); CHKERRQ( ierr );
  ierr = PetscLogStagePush( stage ); CHKERRQ( ierr );

  int Ii, i, j, J;
  PetscScalar v;

  // Set matrix elements for 2D, 5-point stencil in parallel.
  // Each processor needs to insert only elements that it owns locally.
  // Any non-local elements will be sent to the appropriate processor during
  // matrix assembly.
  // Always specify global rows and columns of matrix entries.

  for( Ii = Istart; Ii < Iend; Ii++ ) {
    v = -1.0;
    i = Ii/n;
    j = Ii - i*n;
    //printf("i = %d, j = %d\n", i, j );
    if( i > 0 ) {
      J = Ii - n;
      ierr = MatSetValues( A, 1, &Ii, 1, &J, &v, ADD_VALUES ); CHKERRQ( ierr );
    }
    if( i < m - 1 ) {
      J = Ii + n;
      ierr = MatSetValues( A, 1, &Ii, 1, &J, &v, ADD_VALUES ); CHKERRQ( ierr );
    }
    if( j > 0 ) {
      J = Ii - 1;
      ierr = MatSetValues( A, 1, &Ii, 1, &J, &v, ADD_VALUES ); CHKERRQ( ierr );
    }
    if( j < n - 1 ) {
      J = Ii + 1;
      ierr = MatSetValues( A, 1, &Ii, 1, &J, &v, ADD_VALUES ); CHKERRQ( ierr );
    }
    // Diagonal
    v = 4.0;
    ierr = MatSetValues( A, 1, &Ii, 1, &Ii, &v, ADD_VALUES ); CHKERRQ( ierr );
  }

  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ); CHKERRQ( ierr );
  ierr = PetscLogStagePop(); CHKERRQ( ierr );

  // Set symmetric flag to enable ICC/Cholesky preconditioner
  ierr = MatSetOption( A, MAT_SYMMETRIC, PETSC_TRUE ); CHKERRQ( ierr );

  if( n <= 8 && m <= 8 ) {
    ierr = MatView( A, PETSC_VIEWER_STDOUT_WORLD );
  }

  /*
  PetscViewer viewer;
  PetscDraw draw;
  ierr = PetscViewerCreate( PETSC_COMM_WORLD, &viewer ); CHKERRQ( ierr );
  ierr = PetscViewerSetType( viewer, PETSCVIEWERDRAW ); CHKERRQ( ierr );
  ierr = MatView( A, viewer ); CHKERRQ( ierr );
  ierr = PetscViewerDrawGetDraw( viewer, 0, &draw ); CHKERRQ( ierr );
  ierr = PetscDrawSetPause( draw, 2 ); CHKERRQ( ierr );
  */

  Vec u, b, x;
  //
  ierr = VecCreate( PETSC_COMM_WORLD, &u ); CHKERRQ( ierr );
  ierr = VecSetSizes( u, PETSC_DECIDE, m*n ); CHKERRQ( ierr );
  ierr = VecSetFromOptions(u); CHKERRQ( ierr );
  ierr = VecDuplicate( u, &b ); CHKERRQ( ierr );
  ierr = VecDuplicate( b, &x ); CHKERRQ( ierr );

  PetscBool flg = PETSC_FALSE;
  PetscRandom rctx; // random number generator context
  ierr = PetscOptionsGetBool( NULL, NULL, "-random_exact_sol", &flg, NULL );
  if( flg ) {
    ierr = PetscRandomCreate( PETSC_COMM_WORLD, &rctx ); CHKERRQ( ierr );
    ierr = PetscRandomSetFromOptions( rctx ); CHKERRQ( ierr );
    ierr = VecSetRandom( u, rctx ); CHKERRQ( ierr );
    ierr = PetscRandomDestroy( &rctx ); CHKERRQ( ierr );
  }
  else {
    ierr = VecSet( u, 1.0 ); CHKERRQ( ierr );
  }
  //
  ierr = MatMult( A, u, b ); CHKERRQ( ierr );

  //
  // View exact solution vector if desired
  //
  flg = PETSC_FALSE;
  ierr = PetscOptionsGetBool( NULL, NULL, "-view_exact_sol", &flg, NULL ); CHKERRQ( ierr );
  if( flg ) {
    ierr = VecView( u, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ( ierr );
  }


  // Create linear solver and set various options
  KSP ksp;
  ierr = KSPCreate( PETSC_COMM_WORLD, &ksp ); CHKERRQ( ierr );
  // 
  // set operators
  ierr = KSPSetOperators( ksp, A, A ); CHKERRQ( ierr );
  //
  // various options
  ierr = KSPSetTolerances( ksp, 1.e-2/((m+1)*(n+1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT );
  CHKERRQ( ierr );
  //
  // runtime options:
  // -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  ierr = KSPSetFromOptions( ksp ); CHKERRQ( ierr );


  // Solve the linear system
  ierr = KSPSolve( ksp, b, x ); CHKERRQ( ierr );

  // check error
  // y <- alpha*x + y
  // x <- -u + x
  ierr = VecAXPY( x, -1.0, u ); CHKERRQ( ierr );
  PetscReal norm;
  ierr = VecNorm( x, NORM_2, &norm ); CHKERRQ( ierr );
  PetscInt iters;
  ierr = KSPGetIterationNumber( ksp, &iters ); CHKERRQ( ierr );
  //
  ierr = PetscPrintf( PETSC_COMM_WORLD,
         "Norm of error: %g in %D  iterations\n", (double)norm, iters ); CHKERRQ( ierr );

  ierr = KSPDestroy( &ksp ); CHKERRQ( ierr );
  ierr = VecDestroy( &u ); CHKERRQ( ierr );
  ierr = VecDestroy( &x ); CHKERRQ( ierr );
  ierr = VecDestroy( &b ); CHKERRQ( ierr );
  ierr = MatDestroy( &A ); CHKERRQ( ierr );
  


  PetscPrintf(PETSC_COMM_WORLD, "Program ended normally\n");

  ierr = PetscFinalize(); CHKERRQ( ierr );

  return 0;
}
