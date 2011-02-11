// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/* Computes the gigaflops rate  */
double gflops(int n, double t){
    int64_t n3 = n*n; n3 *= n;
    int64_t flops = (n3*2)/3; 
    double flop_rate = (double) flops / t;
    return ( flop_rate/1.0e9);
}

void Report(int N, int NT, double t0,int matCode,int verif)
{
  double gf = gflops(N,  t0);
  printf("N=%d: %f seconds (%f gflops)\n\n",N,t0,gf);
  printf("#> %d %d %f %f %1c %1c\n\n",N,NT,t0,gf,verif,matCode);
}
