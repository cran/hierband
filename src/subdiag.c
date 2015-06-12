#include <R.h>
#include <Rmath.h>
#include <math.h>
//#include <math.h>

void subdiag_l2norm(double *R, int *dim, double *r) {
  /* 
     Compute the L2 norm of each subdiagonal of a symmetric dim by dim matrix R.

     Assumes R is symmetric.
  */
  int i, j;
  int p = *dim;
  for (i = 0; i < p; i++)
    r[i] = 0;
  for (i = 0; i < p; i++) {
    for (j = 0; j <= i; j++) {
      r[i-j] += R[j + p * i] * R[j + p * i];
    }
  }
  for (i = 0; i < p; i++)
    r[i] = sqrt(r[i]);
}

void formw(int *p, double *w) {
  /*
    Form the w matrix.  For m <= l...
    w[l, m] <- sqrt(2 * l) / (l - m + 1)
   */
  int l, mm, pl;
  double rootl;
  for (l = 0; l < *p - 1; l++) {
    rootl = sqrt( 2 * (l + 1) );
    pl = (*p - 1) * l;
    for (mm = 0; mm <= l; mm++) {
      w[pl + mm] = rootl / (l - mm + 1);
    }
  }
}
