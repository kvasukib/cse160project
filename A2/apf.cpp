/* 
 * Driver for a cardiac elecrophysioly simulatin that uses the
 * Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 */

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <pthread.h>
#include "time.h"
#include "apf.h"
#include "types.h"
using namespace std;


// Utilities
// 

// Allocate a 2D array

_DOUBLE_ **alloc2D(int m,int n){

   _DOUBLE_ **E;
   int nx=n+1, ny=m+1;
   E = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*ny + sizeof(_DOUBLE_)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++) E[j] = (_DOUBLE_*)(E+ny) + j*nx;
   return(E);
}

void init (_DOUBLE_ **E,_DOUBLE_ **E_prev,_DOUBLE_ **R,int m,int n){
    int i,j;
    // Initialization
    for (j=1; j<=m + 1; j++)
        for (i=1; i<= n+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }
    for (j=1; j<=m + 1; j++)
        for (i=n/2+2; i<= n+1 ; i++){
            E_prev[j][i] = 1.0;
    }

    for (j=m/2+2; j<=m+1; j++)
        for (i=1; i<=n+1; i++)
            R[j][i] = 1.0;
}

// External functions
void cmdLine(int argc, char *argv[], _DOUBLE_& T, int& n, int& tx, int& ty, int& do_stats, int& plot_freq);
int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, _DOUBLE_ T, _DOUBLE_ alpha, _DOUBLE_ dt, int do_stats, int plot_freq, int stats_freq);
void printTOD(ofstream& logfile, string mesg);
void ReportStart(ofstream& logfile, _DOUBLE_ dt, _DOUBLE_ T, int m, int n, int tx, int ty);
void ReportEnd(ofstream& logfile, _DOUBLE_ T, int niter, _DOUBLE_ **E_prev, int m,int n, double t0, int tx, int ty);
void *solve_thr( void *arg );

//global vars
pthread_barrier_t barr;
int q, first_q, rest_q; //used to determine partitioning
_DOUBLE_ **E, **R, **E_prev;
 _DOUBLE_ T;
 int m;
 int n;
 int do_stats;
 int plot_freq;
 int tx, ty;

 _DOUBLE_ dt;
 _DOUBLE_ alpha;

int niter;
// This parameter controls the frequncy (in timesteps)
// that summary statistics are reported
// Reduce the value of FREQ to increase the frequency,
// increase the value to raise the frequency
const int STATS_FREQ = 100;

ofstream logfile("Log.txt",ios::out);

int NT = 4;
// Main program
int main(int argc, char** argv)
{
 /*
  *  Solution arrays
  *   E is the "Excitation" variable, a voltage
  *   R is the "Recovery" variable
  *   E_prev is the Excitation variable for the previous timestep,
  *      and is used in time integration
  */


 // Default values for the command line arguments
 T=1500.0;
 m=100;
 n=100;
 do_stats = 0;
 plot_freq = 0;
 tx = 1;
 ty = 1;




// Parse command line arguments
 cmdLine( argc, argv, T, n, tx, ty, do_stats,  plot_freq);
 if (n < 26){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    exit(-1);
 }
 m = n;

 // The log file
 // Do not change the file name or remove this call

   printTOD(logfile, "Simulation begins");

    
 // Allocate contiguous memory for solution arrays
 // The computational box is defined on [1:m+1,1:n+1]
 // We pad the arrays in order to facilitate differencing on the 
 // boundaries of the computation box
 E = alloc2D(m+2,n+2);
 E_prev = alloc2D(m+2,n+2);
 R = alloc2D(m+2,n+2);

 init(E,E_prev,R,m,n);
#ifdef DEBUG
   repNorms(E_prev,-1,dt,m,n,-1,STATS_FREQ);
   if (plot_freq)
     splot(E_prev,-1,-1,m+1,n+1,WAIT);
#endif



 //
 // Initization of various simulation variables
 // Do not change the code these assignments statements, as doing so
 // could cause your submission to be graded incorrectly
 //

 // We compute the timestep dt which determines how long 
 // the code will run for

 // We should always use double precision values for the folowing variables:
 //    rp, dte, dtr, ddt
 //
 // This ensures that the computation of dte and especially dt
 // will not lose precision (i.e. if computed as single precision values)

 if(pthread_barrier_init(&barr, NULL, NT))
 {
    cerr << "could not create barrier\n";
    exit(-1);
 }
 _DOUBLE_ dx = 1.0/n;
 double rp= kk*(b+1)*(b+1)/4;
 double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));void *prime_thr( void *arg );
 double dtr=1/(epsilon+((M1/M2)*rp));
 double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
 dt = (_DOUBLE_) ddt;
 alpha = d*dt/(dx*dx);

 //figure out partitioning
 if(n % NT != 0)
 {
   q = n % NT;
   first_q = (int) ceil(n/NT);
   rest_q = (int) floor(n/NT);

 }
 else
 {
   q = NT;
   first_q = (int) (n+1)/NT;
 }
 // End Initization of various simulation variables

 // Report various information
 // Do not remove this call, it is needed for grading
  ReportStart(logfile, dt, T, m, n, tx, ty);

 // Start the timer
 double t0 = -getTime();

 pthread_t * th_arr = new pthread_t [NT];
 niter = 0;
 for(int i = 0; i < NT; i++)
 {
   int64_t ind = i;
   if(pthread_create(&th_arr[i], NULL, solve_thr, reinterpret_cast<void *> (ind)))
   {
      cerr << "could not create thread " << i;
   }
 }

 cerr << "threads have started\n";
 //int niter = solve(logfile, &E, &E_prev, R, m, n, T, alpha, dt, do_stats, plot_freq,STATS_FREQ);
 for (int t=0; t < NT; t++)
 {
   pthread_join(th_arr[t], NULL);
 }
 t0 += getTime();

 cerr << "about to execute report end\n";
 // Report various information
 // Do not remove this call, it is needed for grading
 ReportEnd(logfile,T,niter,E_prev,m,n,t0,tx,ty);

 cerr << "done with report end\n";

 if (plot_freq){
    printf("\n\nEnter any input to close the program and the plot...");
    int resp;
    scanf("%d",&resp);
  }

 free (E);
 free (E_prev);
 free (R);
}
