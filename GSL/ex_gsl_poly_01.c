#include <stdio.h>
#include <gsl/gsl_poly.h>

int main( int argc, char **argv )
{
  int i;
  double a[6] = {-1, 0, 0, 0, 0, 1};
  double z[10];

  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(6);

  gsl_poly_complex_solve(a, 6, w, z);

  gsl_poly_complex_workspace_free(w);

  // print result
  for( i = 0; i<5; i++ )
  {
    printf("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i+1]);
  }

  return 0;
}

