#include "slepceps.h"

static char help[] =
"Standard symmetric eigenproblem corresponding to the Kinetic\n"
"operator in 1 dimension.\n\n"
"The command line options are:\n"
"  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

int main(int argc, char **argv)
{

  EPS eps;
  EPSType type;
  PetscReal error, tol, re, im;
  PetscScalar kr, ki;
  PetscInt i, Istart, Iend, nev, maxit, its, nconv;

  PetscErrorCode ierr;

  SlepcInitialize( &argc, &argv, (char*)0, help);

  // Set problem size from command line option if provided.
  // The default value is 30
  PetscInt n = 30;
  ierr = PetscOptionsGetInt( NULL, NULL, "-n", &n, NULL ); CHKERRQ(ierr);

  ierr = PetscPrintf( PETSC_COMM_WORLD,
           "\n1-D Laplacian Eigenproblem, n=%D\n\n", n); CHKERRQ(ierr);

  // Create operator matrix that define the eigensystem
  Mat A;
  ierr = MatCreate( PETSC_COMM_WORLD, &A ); CHKERRQ(ierr);
  // Set the size
  ierr = MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, n, n ); CHKERRQ(ierr);
  // Matrix size can be set from command line option, so we need to
  // call these.
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);

  // Start building the matrix
  ierr = MatGetOwnershipRange( A, &Istart, &Iend ); CHKERRQ(ierr);

  for(i = Istart; i < Iend; i++) {
    // lower subdiagonal
    if( i > 0) {
      ierr = MatSetValue( A, i, i-1, 1.0, INSERT_VALUES); CHKERRQ(ierr);
    }
    // upper subdiagonal
    if( i < n-1 ) {
      ierr = MatSetValue( A, i, i+1, 1.0, INSERT_VALUES); CHKERRQ(ierr);
    }
    // diagonal
    ierr = MatSetValue( A, i, i, -2.0, INSERT_VALUES ); CHKERRQ(ierr);
  }
  // Assembly the matrix
  ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
  ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);

  // Eigenvectors (not really used in this case)
  Vec xr, xi;
  ierr = MatCreateVecs( A, NULL, &xr ); CHKERRQ(ierr);
  ierr = MatCreateVecs( A, NULL, &xi ); CHKERRQ(ierr);

  // Create eigensolver and set various options
  ierr = EPSCreate( PETSC_COMM_WORLD, &eps ); CHKERRQ(ierr);

  /* Set operators: standard eigenvalue problem (no overlap matrix) */
  ierr = EPSSetOperators( eps, A, NULL ); CHKERRQ(ierr);

  ierr = EPSSetProblemType( eps, EPS_HEP ); CHKERRQ(ierr);

  // Set solver parameters at runtime
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  // Solve the eigensystem
  ierr = EPSSolve(eps); CHKERRQ(ierr);

  ierr = EPSGetIterationNumber( eps, &its ); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "Number of iterations of the method: %D\n",its); CHKERRQ(ierr);

  ierr = EPSGetType(eps,&type); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type); CHKERRQ(ierr);

  ierr = EPSGetDimensions(eps,&nev,NULL,NULL); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "Number of requested eigenvalues: %D\n", nev); CHKERRQ(ierr);

  ierr = EPSGetTolerances( eps, &tol, &maxit); CHKERRQ(ierr);

  ierr = PetscPrintf( PETSC_COMM_WORLD,
            "Stopping condition: tol=%.4g, maxit=%D\n",
            (double)tol, maxit); CHKERRQ(ierr);

  ierr = EPSGetConverged( eps, &nconv ); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Number of converged eigenpairs: %D\n\n", nconv); CHKERRQ(ierr);

  if ( nconv > 0) {
    // Display eigenvalues and relative errors
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n"); CHKERRQ(ierr);

    for ( i = 0; i < nev; i++) {
      // Get converged eigenpairs: i-th eigenvalue is stored in
      // kr (real part) and ki (imaginary part)
      ierr = EPSGetEigenpair( eps, i, &kr, &ki, xr, xi ); CHKERRQ(ierr);
      // Compute the relative error associated to each eigenpair
      ierr = EPSComputeError( eps, i, EPS_ERROR_RELATIVE, &error); CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
      printf("Pass here ...\n");
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if ( im != 0.0 ) {
        ierr = PetscPrintf( PETSC_COMM_WORLD,
                 " %9f%+9fi %12g\n", (double)re, (double)im, 
                 (double)error); CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf( PETSC_COMM_WORLD,
                  "   %12f       %12g\n", (double)re, (double)error); CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  /*
  Free work space
  */
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&xr); CHKERRQ(ierr);
  ierr = VecDestroy(&xi); CHKERRQ(ierr);
  ierr = SlepcFinalize();

  return ierr;

}
