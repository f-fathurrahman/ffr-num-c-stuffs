
Generate `tags` file for Vim autocompletion:
```
ctags --languages=C -R -f tags ~/WORKS/my_github_repos/ffr-LFDFT-mpi/local/src/petsc-3.8.3
```

## First thing first

Include file (everything):

```c
#include <petsc.h>
```

Skeleton of typical program:

```c
int main( int argc, char** argv ) {
  PetscErrorCode ierr;
  ierr = PetscInitialize( &argc, &argv, NULL, NULL );
  CHKERRQ( ierr );
  // ...
  ierr = PetscFinalize();
  return 0;
}
```


## Working with vectors

Create parallel (MPI) vectors:

```c
Vec x;
PetscInt N = 10;
VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, N, &x );
```

- `VecCreateMPI`
- `VecDestroy`
- `VecDuplicate`
- `VecDuplicateVecs`
- `VecSet`
- `VecSetValues`
- `VecSetRandom`
- `VecView`