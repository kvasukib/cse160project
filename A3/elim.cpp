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

#include <omp.h>
using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

double getTime();

extern double **A, **R;
extern int NT;

void elim(int N)
{
  int i, j, k, Mx;

// If we get stuck, we can count the number of row swaps
// Do this for a small matrix
// int swaps = 0;
  int tid;
  int totalrows = 0;
  int chunk;
#pragma omp parallel private(i,j,k, tid, totalrows,chunk) num_threads(NT)
for ( k = 0; k < N; k++ ) {
//chunk = ceil( (float)(N-k-1)/NT);
//tid = omp_get_thread_num();
//totalrows=0;
/* Partial Pivoting */
#pragma omp single
{
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
}
#pragma omp for// schedule(static, 1)
    for ( i = k+1; i < N; i++ ) 
    {
      A[i][k] /= A[k][k];  
      totalrows++;
      //printf("Thread %d doing %d\n", tid, i);
    }
#pragma omp for// schedule(static, 1)
    for ( i = k+1; i < N; i++ ) {
      double Aik = A[i][k];
      double *Ai = A[i];
      for ( j = k+1; j < N; j++ ) 
        Ai[j] -= Aik * A[k][j];
    }  
 // printf("Thread%d did %d rows\n",tid,totalrows);
  }
//    cout << "Did " << swaps << " row swaps " << endl;
}
