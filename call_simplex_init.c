#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "random.c"

/* Performs initial phase of simplex algorithm to find a feasible point
   Inputs: A, c
   A*x+c <= 0
   s.t. x >= 0 
   Output: x, exitflag (0 = infeasible, 1 = feasible)
*/

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *Ain, *cin, *A, *c, *b, eps, *res, *exitflag, bmax, temp, temp2, xmin;
  int m, n, i, j, k, piv_nb, piv_b, tempi, *basic, *nonbasic, *bpos, nbpos, r, l, found, maxit, it;
  long idum = 1;

  /* read inputs */
  Ain = mxGetPr(prhs[0]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  cin = mxGetPr(prhs[1]);
  eps = mxGetScalar(prhs[2]);
  maxit = mxGetScalar(prhs[3]);
  
  /* output */
  plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
  res = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  exitflag = mxGetPr(plhs[1]);
  
  A = (double*) malloc(sizeof(double)*m*(m+n));
  b = (double*) malloc(sizeof(double)*(m+n));
  c = (double*) malloc(sizeof(double)*m);
  basic = (int*) malloc(sizeof(int)*m);
  nonbasic = (int*) malloc(sizeof(int)*(m+n));

  for (i=0; i<m; i++)
    basic[i] = m+n+i+1;      /* Auxiliary variables (z_i) */
  for (j=0; j<m+n; j++)
    nonbasic[j] = j+1;       /* x_i and slack variables (y_i) */

  /* A contains columns for x_i, y_i */
  for (i=0; i<m; i++)
    if (cin[i] < 0) {
      for (j=0; j<n; j++)
	A[i+m*j] = -Ain[i+m*j];
      for (j=n; j<m+n; j++)
	if (j==n+i)
	  A[i+m*j] = -1;
	else
	  A[i+m*j] = 0;
    }
    else {
      for (j=0; j<n; j++)
	A[i+m*j] = Ain[i+m*j];
      for (j=n; j<m+n; j++)
	if (j==n+i)
	  A[i+m*j] = 1;
	else
	  A[i+m*j] = 0;
    }
  
  for (i=0; i<m; i++)
    if (cin[i] < 0)
      c[i] = -cin[i];
    else
      c[i] = cin[i];
  
  for (j=0; j<n; j++) {       /* Obj. fn = -Sum(Ax+c+y) */
    b[j] = 0;
    for (i=0; i<m; i++)
      b[j] -= A[i+m*j];
  }
  for (j=n; j<m+n; j++)
    if (cin[j-n] < 0)
      b[j] = 1;
    else
      b[j] = -1;

  /* Run simplex method.
     Use random pivot selection. */

  bpos = (int*) malloc(sizeof(int)*(m+n));
  for (it=0; it<maxit; it++) { 

#if 0
    for (i=0;i<m;i++) {
      printf("%.5g \t\t",c[i]);
      for (j=0;j<m+n;j++)
	printf("%.5g \t",A[i+m*j]);
      printf("\n");
    }
    printf("\n");
#endif
#if 0
    for (i=0;i<m+n;i++)
      printf("%.5g \t",b[i]);
    printf("\n\n");
#endif
#if 0
    temp = 0;
    for (i=0;i<m;i++)
      if (basic[i] > m+n)
	temp += c[i];
    printf("obj. = %.5g \n\n",temp);
    if (temp < eps) {
      for (i=0;i<m;i++) {
	printf("%.5g \t\t",c[i]);
	for (j=0;j<m+n;j++)
	  printf("%.5g \t",A[i+m*j]);
	printf("\n");
      }
      printf("\n");
      for (i=0;i<m+n;i++)
	printf("%.5g \t",b[i]);
      printf("\n\n");
    }
#endif

    nbpos = 0;
    for (i=0; i<m+n; i++) 
      if (b[i] > eps && nonbasic[i] <= m+n)	/* ignore exchanged z-variables */      
	for (k=0; k<m; k++)
	  if (A[k+m*i] < -eps) {
	    bpos[nbpos] = i;
	    nbpos++;
	    break;
	  }
    
    /* If the obj. fn. cannot be improved, end the loop */
    if (nbpos == 0)
      break;
    else {
      /* Otherwise, find a random nonbasic pivot variable */
      r = (int)floor((float)nbpos * ran1(&idum));
      piv_nb = bpos[r];
        
      /* Find a basic variable that keeps x_i >= 0 */
      k=0;
      while (k < m) {
	temp = A[k+m*piv_nb];
	if (temp < -eps) {
	  temp2 = -c[k]/temp;
	  xmin = temp2;
	  piv_b = k;
	  k++;
	  break;
	}
	k++;
      }
      while (k < m) {
	temp = A[k+m*piv_nb];
	if (temp < -eps) {
	  temp2 = -c[k]/temp;
	  if (temp2 < xmin) {
	    xmin = temp2;
	    piv_b = k;
	  }
	  else if (temp2 == xmin) 
	    if (ran1(&idum) < 0.5)
	      piv_b = k;
	}
	k++;
      }

#if 0
      /* printf("%.5g \n",A[piv_b+m*piv_nb]); */
      printf("piv %d %d \n\n",basic[piv_b],nonbasic[piv_nb]);
      printf("increase = %.5g \n",-c[piv_b]/A[piv_b+m*piv_nb]*b[piv_nb]);
#endif

      /* Exchange variables */
      temp = -1/A[piv_b+m*piv_nb];
        
      for (j=0; j<piv_nb; j++) { 
	A[piv_b+m*j] *= temp;
	if (fabs(A[piv_b+m*j]) < eps)
	  A[piv_b+m*j] = 0;
      }
      for (j=piv_nb+1; j<m+n; j++) {
	A[piv_b+m*j] *= temp;
	if (fabs(A[piv_b+m*j]) < eps)
	  A[piv_b+m*j] = 0;
      }
      A[piv_b+m*piv_nb] = -temp;
      if (fabs(A[piv_b+m*piv_nb]) < eps)
	A[piv_b+m*piv_nb] = 0;
      c[piv_b] *= temp;
      if (fabs(c[piv_b]) < eps)
	c[piv_b] = 0;
      
      for (i=0; i<piv_b; i++) {
	temp = A[i+m*piv_nb];
	for (j=0; j<m+n; j++) {
	  A[i+m*j] += temp*A[piv_b+m*j];
	  if (fabs(A[i+m*j]) < eps)
	    A[i+m*j] = 0;
	}
	c[i] += temp*c[piv_b];
	if (fabs(c[i]) < eps)
	  c[i] = 0;
	A[i+m*piv_nb] -= temp;
	if (fabs(A[i+m*piv_nb]) < eps)
	  A[i+m*piv_nb] = 0;
      }
      for (i=piv_b+1; i<m; i++) {
	temp = A[i+m*piv_nb];
	for (j=0; j<m+n; j++) {
	  A[i+m*j] += temp*A[piv_b+m*j];
	  if (fabs(A[i+m*j]) < eps)
	    A[i+m*j] = 0;
	}
	c[i] += temp*c[piv_b];
	if (fabs(c[i]) < eps)
	  c[i] = 0;
	A[i+m*piv_nb] -= temp;
	if (fabs(A[i+m*piv_nb]) < eps)
	  A[i+m*piv_nb] = 0;
      }
      
      temp = b[piv_nb];
      for (j=0; j<m+n; j++) {
	b[j] += temp*A[piv_b+m*j];
	if (fabs(b[j]) < eps)
	  b[j] = 0;
      }
      b[piv_nb] -= temp;
      if (fabs(b[piv_nb]) < eps)
	b[piv_nb] = 0;
      
      tempi = basic[piv_b];
      basic[piv_b] = nonbasic[piv_nb];
      nonbasic[piv_nb] = tempi;
    }
  }
  
  *exitflag = 1;
  for (i=0; i<m; i++)
    if (basic[i] > m+n && c[i] > eps) {
      *exitflag = 0;
      break;
    }

  if (*exitflag) {
    /* Exchange remaining nonbasic x_i with z_i */
    for (i=0; i<m+n; i++)
      if (nonbasic[i] <= n) {
    	found = 0;
        for (k=0; k<m; k++) 
	  if (basic[k] > m+n && A[k+m*i] < -eps) {
	    piv_b = k;
	    piv_nb = i;
	    /* Exchange variables */
	    temp = -1/A[piv_b+m*piv_nb];
	    
	    for (j=0; j<piv_nb; j++)
	      A[piv_b+m*j] *= temp;
	    for (j=piv_nb+1; j<m+n; j++)
	      A[piv_b+m*j] *= temp;
	    A[piv_b+m*piv_nb] = -temp;
	    c[piv_b] *= temp;
	    
	    for (i=0; i<piv_b; i++) {
	      temp = A[i+m*piv_nb];
	      for (j=0; j<m+n; j++)
                A[i+m*j] += temp*A[piv_b+m*j];
	      c[i] += temp*c[piv_b];
	      A[i+m*piv_nb] -= temp;
	    }
	    for (i=piv_b+1; i<m; i++) {
	      temp = A[i+m*piv_nb];
	      for (j=0; j<m+n; j++)
                A[i+m*j] += temp*A[piv_b+m*j];
	      c[i] += temp*c[piv_b];
	      A[i+m*piv_nb] -= temp;
	    }
        
	    temp = b[piv_nb];
	    for (j=0; j<m+n; j++)
	      b[j] += temp*A[piv_b+m*j];
	    b[piv_nb] -= temp;
	    
	    tempi = basic[piv_b];
	    basic[piv_b] = nonbasic[piv_nb];
	    nonbasic[piv_nb] = tempi;

	    found = 1;
	    break;
	  }
	if (!found)
	  res[nonbasic[i]-1] = 0;
      }
  }
    
  for (i=0; i<m; i++)
    if (basic[i] <= n)
      res[basic[i]-1] = c[i];
  
  free(A);
  free(b);
  free(c);
  free(basic);
  free(nonbasic);
  free(bpos); 
}
