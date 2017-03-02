#include "slepcsys.h"

static char help[] = "Simple Hello World example program in SLEPC\n";

int main( int argc, char **argv)
{
  int ierr;

  SlepcInitialize( &argc, &argv, (char*)0, help);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Hello World SLEPC by ffr\n");
  CHKERRQ(ierr);

  ierr = SlepcFinalize();

  return ierr;
}
