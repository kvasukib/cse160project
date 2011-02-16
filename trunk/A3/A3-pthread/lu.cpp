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
#include <pthread.h>
#include "time.h"

using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

double getTime();
int verify(int N);
pthread_barrier_t barr;

double **A, **R;

void Report(int N, int NT, double t0,int matCode,int verif);

void initTimer(int NT);
void FinalizeTimer();
void printMat(int N,double **A){
  int i,j ;
  if (N < 10){
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%8.3g ",A[i][j]);
        }
        printf("\n"); printf("\n");
    }
    printf("\n");
 }
}
int NT;
int N;
double t0;
// External routines
void * elim_thr(void * arg);
void initMat(int N, double **A, double **R, int matCode);
void cmdLine(int argc, char *argv[], int& n, int& nt, int& matCode);

int main(int argc,char *argv[])
{
  int i, matCode;

  cmdLine(argc, argv, N, NT, matCode);
  printf("N = %d, nt = %d\n",N,NT);

  switch (matCode){
      case 'L':
          cout << "Using the '5-point' Matrix" << endl;
          break;

      case 'H':
          cout << "Using a 'Hadamard' Matrix" << endl;
          break;

      case 'D':
          cout << "Using a simple matrix 2*I" << endl;
          break;

      default:
          cerr << "No such matrix." << endl;
          break;
  }

  /* Allocate A and  R as N by N matrices */
  assert(A = (double **)malloc(N*sizeof(double *)));
  assert(R = (double **)malloc(N*sizeof(double *)));
  for (i = 0; i < N; i++) {
    assert(A[i] = (double *)malloc(N*sizeof(double)));
    assert(R[i] = (double *)malloc(N*sizeof(double)));
  } 

printf("# of pthreads: %d\n",NT);

 if(pthread_barrier_init(&barr, NULL, NT))
 {
    cerr << "could not create barrier\n";
    exit(-1);
 }


  initMat(N,A, R,matCode);
  printMat(N,A);

//warmup

  pthread_t * th_arr_warm = new pthread_t [NT];
initTimer(NT);
for(int i = 0; i < NT; i++)
 {
   int64_t ind = i;
   if(pthread_create(&th_arr_warm[i], NULL, elim_thr, reinterpret_cast<void *> (ind)))
   {
      cerr << "could not create thread " << i;
   }
 }

 for (int t=0; t < NT; t++)
 {
   pthread_join(th_arr_warm[t], NULL);
 }
 FinalizeTimer();

  initMat(N,A, R,matCode);
  pthread_t * th_arr = new pthread_t [NT];

 initTimer(NT);
 for(int i = 0; i < NT; i++)
 {
   int64_t ind = i;
   if(pthread_create(&th_arr[i], NULL, elim_thr, reinterpret_cast<void *> (ind)))
   {
      cerr << "could not create thread " << i;
   }
 }

 for (int t=0; t < NT; t++)
 {
   pthread_join(th_arr[t], NULL);
 }
 FinalizeTimer();

  /* print results */
  printMat(N,A);
//  printf("N=%d: %f seconds\n",N,t0);
//  double t1 = -getTime();
  int verif = verify(N);
//  t1 += getTime();
//  cout << "Verification time: " << t1 << " seconds\n";

    Report(N, NT,t0,matCode,verif);
  return(0);
}

/* EOF ludcmp.c */
