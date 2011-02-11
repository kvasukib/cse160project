// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
//
// Scott B. Baden, UCSD, 2/2/11
//

#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <math.h>

#ifdef OPENMP_
#include <omp.h>
#endif
using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

extern double **A, **R;

int verify(int N){
  /* Verify */
  int i, j, k;
  double epsilon = 1.0e-6;
  double *b = (double *) malloc(sizeof(double)*N); assert(b);
  double *x = (double *) malloc(sizeof(double)*N); assert(x);
  double *y = (double *) malloc(sizeof(double)*N); assert(y);
  double *save = (double *) malloc(sizeof(double)*N); assert(save);
  // Choose a simple RHS
  for (i = 0; i < N; i++)
      b[i] = 1.0;

  // Verify ...
  // We are solving Ax = b ->  LUx = b 
  //   First we solve Ly = b, then Ux = y
  // We solve the first eqn with forward subtitution,
  // the second with backward subsitution
  // When done, verify that the solution is correct by computing 
  // the residual  r = Ax - b
  // 

  // Our storage scheme overwrites the diagonal values of L
  // with the diagonal values of U
  // Fortunately, the diagonal values of L are all ones,
  // so we make a copy of the U diagonal and restore it when we need it
  
    for (i = 0; i < N; i++){
        save[i] = A[i][i];
        A[i][i] = 1.0;
    }
  // Forward Solve
    y[0] = b[0] / A[0][0];
    for (i = 1; i < N; i++){
        double sum = 0.0;
        for (k = 0; k < i; k++)
            sum += A[i][k]*y[k];
        y[i] = (b[i] - sum)/A[i][i];
    }

// We've now solved Ly = b
// Next, we solve Ux = y, but first we have to restore
// the diagonal values we saved away 
    for (i = 0; i < N; i++){
        A[i][i]  = save[i];
    }

    x[N-1] = y[N-1]/A[N-1][N-1];
    for (i = N-2; i >= 0; i--) {
        double sum = 0.0;
        for (k = i+1; k < N; k++)
            sum += A[i][k]*x[k];
        x[i] = (y[i] - sum)/A[i][i];
    }

//    Compute r = Ax-b
    double r = 0;
    for (i=0; i<N; i++){
       double yi = 0.0;
       for (j=0; j<N; j++)
          yi += R[i][j] * x[j];
       double ri = yi - b[i];
       r += ri*ri;
    }
    double resid = sqrt(r)/N;
    if (resid < epsilon){
        cout << "LU decomposition verified OK [eps < " << epsilon << "]\n";
        return ('Y');
    }
    cout << "Non-zero residual = " << resid << endl;
    return ('N');
}
