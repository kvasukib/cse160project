// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//

//
// Initialize the matrix.
// Scott B. Baden, UCSD, 2/2/11
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#ifdef OPENMP_
#include <omp.h>
#endif
using namespace std;

void initLMat(int N, double **A, double **R);
void initRMat(int N, double **A, double **R);
void initHMat(int N, double **A, double **R);

void initMat(int N, double **A, double **R, int matCode){

  switch (matCode){
      case 'L':              // '5-point' Matrix, does not need partial pivoting
          initLMat(N, A, R);
          break;

      case 'H':
          initHMat(N, A, R); // 'Hadamard' Matrix, requires pivoting
          break;

      case 'D':              // Simple matrix 2*I, I = Identity Matrix
          initRMat(N, A, R);
          break;

      default:
          cerr << "No such matrix." << endl;
          break;
  }
}
void initLMat(int N, double **A, double **R){
// Generates an n^2-by-n^2  matrix corresponding
// to the discretization of the 5-point Laplacian
// on a n x n mesh; A is symmetric positive definite 


 int isqrt = sqrt(N);
 if (isqrt*isqrt != N){
     cerr << "\n *** N is not an integral square root. Exiting" << endl << endl;
     exit(-1);
 }
 int n2 = N;
 int n = isqrt;
 int i, j;
 for (i=0; i<n2; i++)
    for (j=0;j<n2;j++)
        A[i][j] = 0;


 for (i=0; i<n2; i++)
     A[i][i] = 4;
 for (i=0; i<n; i++)
    for (j=0;j<n-1;j++){
       int x = n*i+j;
       A[x][x+1]=-1;
       A[x+1][x]=-1;
    }
 for (i=n; i<n2; i++){
    A[i][i-n]=-1;
    A[i-n][i]=-1;
 }

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      R[i][j] = A[i][j];

  return;
}

// A simpler input matrix

void initRMat(int N, double **A, double **R){
  int i,j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = 0.0;
      if (i == j) A[i][j] = 2.0;
      R[i][j] = A[i][j];
    }
  }
  return;
}

// 
// Generate a dyadic(Paley) ordered Hadamard matrix
// N must be a power of 2.
// Obtained from
// www.mathworks.com/matlabcentral/fileexchange/7158-dyadic-paley-ordered-hadamard-matrix
// This matrix is different from that returned by the Matlab "hadamard" function
//

void mtxdya(int N, double **D){

    int q = log2((double) N);
    int N2 = 2 << (q-1);
    if (N2 != N){
        cerr << "N must be an integer power of 2" << endl;
        exit(-1);
    }


    for (int u=1; u<N+1; u++){           // for u = 1:N
        int binu = u-1;        // dec2bin(u-1,q);
        for (int v=1; v<N+1; v++){       // for v = 1:N
            int binv = v-1;       // dec2bin(v-1,q);
            int temp = 0;
            for (int i=1; i<q+1; i++){   // for i = 1:q
                // bin2dec(binu(i))*bin2dec(binv(q+1-i));
                temp += (((binu >> (q-i)) & 1) & ((binv >> (i-1)) & 1) );
            }

            int s = temp & 1; // D[u][v]=(-1)^temp;
            D[u-1][v-1]  = (double) (1 - (s<<1));
        }
    }
}

void initHMat(int N, double **A, double **R){
    mtxdya(N, A);
    int i,j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            R[i][j] = A[i][j];
        }
  }
  return;
}
