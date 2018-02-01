#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cusparse.h>
#include <cublas_v2.h>


void genLaplace( int *row_ptr, int *col_ind, float *val,
                 int M, int N, int nz, float *rhs )
{
  if( M != N ) {
    printf("ERROR M != N : %d != %d\n", M, N );
    exit(EXIT_FAILURE);
  }

  int n = (int)sqrt((double)N);

  printf("Laplace dimension = %d\n", n);

  int idx = 0;
  for( int i = 0; i < N; i++ ) {
    int ix = i%n;
    int iy = i%n;
    row_ptr[i] = idx;
    // up
    if( iy > 0 ) {
      val[idx] = 1.0;
      col_ind[idx] = i - n;
      idx++;
    }
    else {
      rhs[i] -= 1.0;
    }
    // left
    if( ix > 0 ) {
      val[idx] = 1.0;
      col_ind[idx] = i - 1;
      idx++;
    }
    else {
      rhs[i] -= 0.0;
    }
    // center
    val[idx] = -4.0;
    col_ind[idx] = i;
    idx++;
    //right
    if( ix < n-1 ) {
      val[idx] = 1.0;
      col_ind[idx] = i+1;
      idx++;
    }
    else {
      rhs[i] -= 0.0;
    }
    //down
    if( iy < n-1) {
      val[idx] = 1.0;
      col_ind[idx] = i + n;
      idx++;
    }
    else {
      rhs[i] -= 0.0;
    }
  }
  row_ptr[N] = idx;
}

int main( int argc, char** argv )
{
  const int MAX_ITER = 1000;
  
  printf("conjGrad starting ...\n");

  int qatest = 0;

  cudaDeviceProp deviceProp;
  //int devID = findCudaDevice( argc, (const char**)argv );
  int devID = 0;
  printf("GPU selected Device ID = %d\n", devID);

  cudaGetDeviceProperties( &deviceProp, devID );
  printf("GPU device has %d multiprocessors, SM %d.%d compute capabilities\n\n",
         deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);

  int M, N;
  int nz;
  int *I, *J;
  float *val, *x, *rhs;

  // Generate a random tridiagonal symmetric matrix in CSR format
  M = N = 16384;
  nz = 5*N-4*(int)sqrt((double)N);
  I = (int*)malloc(sizeof(int)*(N+1)); // CSR row pointers for matrix A
  J = (int*)malloc(sizeof(int)*nz);  // CSR column indices for matrix A
  val = (float*)malloc(sizeof(float)*nz); // CSR values for matrix A
  x = (float*)malloc(sizeof(float)*N);
  rhs = (float*)malloc(sizeof(float)*N);

  // initialize RHS
  int i;
  for( i = 0; i < N; i++ ) {
    rhs[i] = 0.0;
    x[i] = 0.0;
  }

  genLaplace( I, J, val, M, N, nz, rhs );

  // CUBLAS context
  cublasHandle_t cublasHandle = 0;
  cublasStatus_t cublasStatus;
  cublasStatus = cublasCreate( &cublasHandle );
//  printf("Creating CUBLAS context: cublasStatus = %d\n", cublasStatus);

  // XXX: should check cublasStatus

  // CUSPARSE context
  cusparseHandle_t cusparseHandle = 0;
  cusparseStatus_t cusparseStatus;
  cusparseStatus = cusparseCreate( &cusparseHandle );

  // XXX: should check cusparseStatus

  // description of matrix A
  cusparseMatDescr_t descr = 0;
  cusparseStatus = cusparseCreateMatDescr( &descr );

  // define properties of the matrix
  cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL );
  cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO );

  int *d_col, *d_row;
  float *d_val;
  float *d_x, *d_y, *d_r, *d_p, *d_omega;

  // allocate device memory
  cudaMalloc( (void**)&d_col, nz*sizeof(int) );
  cudaMalloc( (void**)&d_row, (N+1)*sizeof(int) );
  cudaMalloc( (void**)&d_val, nz*sizeof(float) );
  cudaMalloc( (void**)&d_x, N*sizeof(float) );
  cudaMalloc( (void**)&d_y, N*sizeof(float) );
  cudaMalloc( (void**)&d_r, N*sizeof(float) );
  cudaMalloc( (void**)&d_p, N*sizeof(float) );
  cudaMalloc( (void**)&d_omega, N*sizeof(float) );


  printf("\nconjGrad ends normally.\n");
}

