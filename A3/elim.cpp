/*
* Unblocked LU decomposition
* based on a parallel code written by
*  Anu Rao 11/1/94 CSC 652
*
* Code modifications made by 
* Scott B. Baden, UCSD, 1/27/11
*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef OPENMP_
#include <omp.h>
#endif
using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

double getTime();

extern double **A, **R;


void elim(int N)
{
  int i, j, k, Mx;

// If we get stuck, we can count the number of row swaps
// Do this for a small matrix
// int swaps = 0;

  for ( k = 0; k < N; k++ ) {
/* Partial Pivoting */
      Mx = k;
      for ( i = k+1; i < N; i++ ) {
	  if (fabs(A[i][k]) > fabs(A[Mx][k]))
	     Mx = i;
      }
      if (Mx > k){
//        swaps++;
         for ( j = k; j < N; j++ ){
             double t = A[Mx][j];
             A[Mx][j] = A[k][j];
             A[k][j] = t;
         }
      }
/* End Partial Pivoting */

    for ( i = k+1; i < N; i++ ) 
      A[i][k] /= A[k][k];  

    for ( i = k+1; i < N; i++ ) {
      double Aik = A[i][k];
      double *Ai = A[i];
      for ( j = k+1; j < N; j++ ) 
        Ai[j] -= Aik * A[k][j];
    }  
  }
//    cout << "Did " << swaps << " row swaps " << endl;
}
