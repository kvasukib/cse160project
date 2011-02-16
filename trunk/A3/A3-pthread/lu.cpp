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
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#define max(a,b) ((a) > (b) ? (a) : (b))

double getTime();
int verify(int N);
pthread_barrier_t barr;

double **A, **R;

void Report(int N, int NT, double t0,int matCode,int verif);


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
// External routines
void * elim_thr(void * arg);
void initMat(int N, double **A, double **R, int matCode);
void cmdLine(int argc, char *argv[], int& n, int& nt, int& matCode);

int main(int argc,char *argv[])
{
  int i, matCode;
  double t0;

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

#ifdef _OPENMP
    omp_set_num_threads(NT);
#pragma omp parallel
#pragma omp single
    printf("# of openMP threads: %d\n",NT);
#else
    printf("Serial run (no openmp)\n"); 
#endif

 if(pthread_barrier_init(&barr, NULL, NT))
 {
    cerr << "could not create barrier\n";
    exit(-1);
 }


  initMat(N,A, R,matCode);
  printMat(N,A);

  pthread_t * th_arr = new pthread_t [NT];

  t0 = -getTime();
 for(int i = 0; i < NT; i++)
 {
   int64_t ind = i;
   if(pthread_create(&th_arr[i], NULL, elim_thr, reinterpret_cast<void *> (ind)))
   {
      cerr << "could not create thread " << i;
   }
 }

  //Warm up the code by running LU once, before we collect timings.
  //elim(N);
 //
  // Don't forget to re-initialize A and R */
  //initMat(N, A, R, matCode);

  //elim(N);

 for (int t=0; t < NT; t++)
 {
   pthread_join(th_arr[t], NULL);
 }
  t0 += getTime();

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
