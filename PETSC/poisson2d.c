static char help[] = "Solves 2D Poisson equation using multigrid.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>

extern PetscErrorCode ComputeJacobian(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);
extern PetscErrorCode ComputeTrueSolution(DM, Vec);
extern PetscErrorCode VecView_VTK(Vec, const char [], const char []);

typedef enum {DIRICHLET, NEUMANN} BCType;

typedef struct {
  PetscScalar uu, tt;
  BCType      bcType;
} UserContext;

int main(int argc,char **argv)
{
  KSP            ksp;
  DM             da;
  UserContext    user;
  PetscInt       bc;

  PetscInitialize(&argc,&argv,(char*)0,help);
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,-11,-11,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);

  KSPSetDM(ksp,(DM)da);
  DMSetApplicationContext(da,&user);

  user.uu     = 1.0;
  user.tt     = 1.0;
  bc          = (PetscInt)NEUMANN; /* Use Neumann Boundary Conditions */
  user.bcType = (BCType)bc;


  KSPSetComputeRHS(ksp,ComputeRHS,&user);
  KSPSetComputeOperators(ksp,ComputeJacobian,&user);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,NULL,NULL);

  DMDestroy(&da);
  KSPDestroy(&ksp);
  PetscFinalize();
  return 0;
}

PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscInt       i,j,M,N,xm,ym,xs,ys;
  PetscScalar    Hx,Hy,pi,uu,tt;
  PetscScalar    **array;
  DM             da;

  KSPGetDM(ksp,&da);
  DMDAGetInfo(da, 0, &M, &N, 0,0,0,0,0,0,0,0,0,0);
  uu   = user->uu; tt = user->tt;
  pi   = 4*atan(1.0);
  Hx   = 1.0/(PetscReal)(M);
  Hy   = 1.0/(PetscReal)(N);

  DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0); /* Fine grid */
  /* printf(" M N: %d %d; xm ym: %d %d; xs ys: %d %d\n",M,N,xm,ym,xs,ys); */
  DMDAVecGetArray(da, b, &array);
  for (j=ys; j<ys+ym; j++) {
   for (i=xs; i<xs+xm; i++) {
     array[j][i] = -PetscCosScalar(uu*pi*((PetscReal)i+0.5)*Hx)*PetscCosScalar(tt*pi*((PetscReal)j+0.5)*Hy)*Hx*Hy;
  }
}
DMDAVecRestoreArray(da, b, &array);
VecAssemblyBegin(b);
VecAssemblyEnd(b);

/* force right hand side to be consistent for singular matrix */
/* note this is really a hack, normally the model would provide you with a consistent right handside */
if (user->bcType == NEUMANN) {
  MatNullSpace nullspace;

  MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
  MatNullSpaceRemove(nullspace,b);
  MatNullSpaceDestroy(&nullspace);
}
return(0);
}

PetscErrorCode ComputeJacobian(KSP ksp,Mat J, Mat jac,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscInt       i, j, M, N, xm, ym, xs, ys, num, numi, numj;
  PetscScalar    v[5], Hx, Hy, HydHx, HxdHy;
  MatStencil     row, col[5];
  DM             da;

  KSPGetDM(ksp,&da);
  DMDAGetInfo(da,0,&M,&N,0,0,0,0,0,0,0,0,0,0);
  Hx    = 1.0 / (PetscReal)(M);
  Hy    = 1.0 / (PetscReal)(N);
  HxdHy = Hx/Hy;
  HydHx = Hy/Hx;
  DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);
  for (j=ys; j<ys+ym; j++) {
   for (i=xs; i<xs+xm; i++) {
     row.i = i; row.j = j;

   if (i==0 || j==0 || i==M-1 || j==N-1) {
     if (user->bcType == DIRICHLET) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Dirichlet boundary conditions not supported !\n");
 else if (user->bcType == NEUMANN) {
   num=0; numi=0; numj=0;
   if (j!=0) {
     v[num] = -HxdHy;              col[num].i = i;   col[num].j = j-1;
     num++; numj++;
   }
   if (i!=0) {
     v[num] = -HydHx;              col[num].i = i-1; col[num].j = j;
     num++; numi++;
   }
   if (i!=M-1) {
     v[num] = -HydHx;              col[num].i = i+1; col[num].j = j;
     num++; numi++;
   }
   if (j!=N-1) {
     v[num] = -HxdHy;              col[num].i = i;   col[num].j = j+1;
     num++; numj++;
   }
   v[num] = ((PetscReal)(numj)*HxdHy + (PetscReal)(numi)*HydHx); col[num].i = i;   col[num].j = j;
         num++;
         MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);
       }
     } else {
       v[0] = -HxdHy;              col[0].i = i;   col[0].j = j-1;
       v[1] = -HydHx;              col[1].i = i-1; col[1].j = j;
       v[2] = 2.0*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
       v[3] = -HydHx;              col[3].i = i+1; col[3].j = j;
       v[4] = -HxdHy;              col[4].i = i;   col[4].j = j+1;
       MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);
     }
   }
 }
 MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
 MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
 if (user->bcType == NEUMANN) {
   MatNullSpace nullspace;

   MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
   MatSetNullSpace(J,nullspace);
   MatNullSpaceDestroy(&nullspace);
 }
 return(0);
}
