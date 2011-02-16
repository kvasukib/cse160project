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

using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

double getTime();
extern pthread_barrier_t barr;
extern double **A, **R;
extern int NT;
extern int N;

int q;
int first_q;
int rest_q;

int totalwork;

void * elim_thr(void * arg)
{
  int i, j, k, Mx;
  int64_t id = reinterpret_cast<int64_t>(arg);
  int tid = id;
// If we get stuck, we can count the number of row swaps
// Do this for a small matrix
// int swaps = 0;
for ( k = 0; k < N; k++ ) {

  if(tid==0)
  {
    totalwork = N - (k+1);
    if(totalwork % NT != 0)
    {
      q = totalwork % NT;
      first_q = ceil( (float) totalwork/NT);
      rest_q = floor( (float) totalwork/NT);
    }
    else
    {
      q = NT;
      first_q = (int) totalwork/NT;
      
    }
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
    }
/* End Partial Pivoting */

   int rc = pthread_barrier_wait(&barr);
   if (rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
   {
      cerr << "could not wait on barrier\n";
      exit(-1);
   }

  int start;
  int end;
  if(tid < q)
  {
    start = k + 1 + id * first_q;
    end = start + first_q;
  }
  else
  {
    start = k + 1 + (q*first_q) + (id-q)*rest_q;
    end = start+ rest_q;
  }  
    
    for ( i = start; i < end; i++ ) 
    {
      A[i][k] /= A[k][k];  
      //totalrows++;
      //printf("Thread %d doing %d\n", tid, i);
    }

   rc = pthread_barrier_wait(&barr);
   if (rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
   {
      cerr << "could not wait on barrier\n";
      exit(-1);
   }

    for ( i = start; i < end; i++ ) {
      double Aik = A[i][k];
      double *Ai = A[i];
      for ( j = k+1; j < N; j++ ) 
        Ai[j] -= Aik * A[k][j];
    }  
 // printf("Thread%d did %d rows\n",tid,totalrows);
//    cout << "Did " << swaps << " row swaps " << endl;
}
  pthread_exit(NULL);
  return 0;
}
