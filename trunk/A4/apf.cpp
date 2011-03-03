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
#include "time.h"
#include "apf.h"
#include "types.h"
#include "mpi.h"
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

void init (_DOUBLE_ **E,_DOUBLE_ **E_prev,_DOUBLE_ **R,int size,int m,int n, int rank, int tx, int ty){
    int i,j;


 for (j=1; j<=m+1; j++)
   for (i=1; i<=n+1; i++)
     E_prev[j][i] = R[j][i] = 0.0;


 if (size == 1){//single processor
     for (j=1;j<=m+1;j++)
         for (i=n/2+2;i<=n+1;i++)
             E_prev[j][i] = 1.0;
     for (j=m/2+2;j<=m+1;j++)
         for (i=1;i<=n+1;i++)
             R[j][i] = 1.0;
 }
 else if (tx == 1){//horizontal slices
     for (j=1; j<=m+1; j++)
         for (i=n/2+2; i<=n+1; i++)
             E_prev[j][i] = 1.0;
     if (rank == size/2-1)
         for (i = 1; i <= n+1;i++)
             R[m+2][i] = 1.0;
     else if (rank == size/2){//top half
         for (j=1;j <= m+2;j++)
             for (i = 1; i <= n+1;i++)
                 R[j][i] = 1.0;
     } else if (rank > size/2)
         for (j = 0;j<=m+2;j++)
             for (i = 1; i <= n+1;i++)
                 R[j][i] = 1.0;
 } else if (ty == 1){//Split up vertically
    for (j = m/2+2;j <= m+1;j++)
        for (i = 0;i <= n+2;i++)
            R[j][i] = 1.0;
    if (rank == 0)
        for (j = m/2+2;j <= m+1;j++)
            R[j][0] = 0.0;
    if (rank == size - 1)
        for (j = m/2+2;j <= m+1;j++)
            R[j][n+2] = 0.0;
    if (rank == size/2-1)
        for (j = 1; j <= m+1;j++)
            E_prev[j][n+2] = 1.0;
    if (rank == size/2){//top half
        for (j = 1;j <= m+1;j++)
            for (i = 1;i <= n+2;i++)
                E_prev[j][i] = 1.0;
    }
    if (rank > size/2){
        for (j = 1;j <= m+1;j++)
            for (i = 0;i <= n+2;i++)
                E_prev[j][i] = 1.0;
     }
 } else {
      printf("No 2 dimensional support");
      return;
 }



/*
    // Initialization
    for (j=1; j<=my_m + 1; j++)
        for (i=1; i<= my_n+1; i++){
            E_prev[j][i] = R[j][i] = 0;
    }
    for (j=1; j<=my_m + 1; j++){
      // for (i=n/2+2; i<= n+1 ; i++){
       for (i=n/2+2; i<=my_n+1; i++){
            //if(rank*my_n + i >= n/2+2)
              E_prev[j][i] = 1.0;
       }
    }
   
   // for (j=m/2+2; j<=m+1; j++){
    for (j=1; j <= my_m+1; j++){
        if(rank*my_m + j >= n/2+2){
            for (i=1; i<=my_n+1; i++)
                R[j][i] = 1.0;
        }
    }
if(rank==1){
    for (j=1; j <= my_m+1; j++){
            for (i=1; i<=my_n+1; i++)
                cerr << R[j][i];
        cerr << '\n';
        }
    }

char c;
cin >> c;*/
}

// External functions
void cmdLine(int argc, char *argv[], _DOUBLE_& T, int& n, int& tx, int& ty, int& do_stats, int& plot_freq, int& noComm);
int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, _DOUBLE_ T, _DOUBLE_ alpha, _DOUBLE_ dt, int do_stats, int plot_freq, int stats_freq, _DOUBLE_ *** tmp_entire, int rank, int size, int full_n, int tx, int ty, int noComm);
void printTOD(ofstream& logfile, string mesg);
void ReportStart(ofstream& logfile, _DOUBLE_ dt, _DOUBLE_ T, int m, int n, int tx, int ty, int noComm);
void ReportEnd(ofstream& logfile, _DOUBLE_ T, int niter, _DOUBLE_ **E_prev, int m,int n, double t0, int tx, int ty);
// Main program
int main(int argc, char** argv)
{
  MPI_Init( &argc, &argv);
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Process %d of %d reporting\n", rank, size);
 /*
  *  Solution arrays
  *   E is the "Excitation" variable, a voltage
  *   R is the "Recovery" variable
  *   E_prev is the Excitation variable for the previous timestep,
  *      and is used in time integration
  */
 _DOUBLE_ **E, **R, **E_prev;
 _DOUBLE_ ** tmp_entire;

 // Default values for the command line arguments
 _DOUBLE_ T=1500.0;
 int m=100,n=100;
 int do_stats = 0;
 int plot_freq = 0;
 int tx = 1, ty = 1;
 int noComm = 0;
// This parameter controls the frequncy (in timesteps)
// that summary statistics are reported
// Reduce the value of FREQ to increase the frequency,
// increase the value to raise the frequency
 const int STATS_FREQ = 1;

// Parse command line arguments
 cmdLine( argc, argv, T, n, tx, ty, do_stats,  plot_freq, noComm);
 if (n < 26){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    exit(-1);
 }
 m = n;

 // The log file
 // Do not change the file name or remove this call
   ofstream logfile("Log.txt",ios::out);
   if(rank==0){
   printTOD(logfile, "Simulation begins");
   }

  //calculate number of rows for each proc
   int rows_per_thread;
   int my_m;
   int my_n;

   if((n+1) % size == 0){
     rows_per_thread = (int) ceil((n+1)/size);
     my_m = (n+1)/ty -1;
     my_n = (n+1)/tx -1;
   }
   else
   { 
     return 1;
   }

  if(size ==1)
    noComm =1;
 // Allocate contiguous memory for solution arrays
 // The computational box is defined on [1:m+1,1:n+1]
 // We pad the arrays in order to facilitate differencing on the 
 // boundaries of the computation box
   E = alloc2D(my_m+2,my_n+2);
   E_prev = alloc2D(my_m+2,my_n+2);
   R = alloc2D(my_m+2,my_n+2);

 tmp_entire = alloc2D(m+2, n+2);
 init(E,E_prev,R,size,my_m,my_n,rank,tx,ty);
 /*if(1)
 {
   tmp_entire = alloc2D(m+2, n+2);
  
    int i,j;
    // Initialization
    for (j=1; j<=m + 1; j++)
        for (i=1; i<= n+1; i++){
            tmp_entire[j][i] = 0;
    }
    for (j=1; j<=m + 1; j++){
       for (i=n/2+2; i<= n+1 ; i++){
              tmp_entire[j][i] = 1.0;
       }
    }*/
/*
     for(i = 1; i <=n+1; i++)
     {
        for(j=1; j <=n+1; j++)
            printf("Reading i=%d j=%d\n", i, j);
            cerr << tmp_entire[i][j];
        cerr << '\n';
     }
*/
#ifdef DEBUG
 if(rank==0){
   repNorms(tmp_entire,-1,dt,m,n,-1,STATS_FREQ);
   if (plot_freq)
     splot(tmp_entire,-1,-1,m+1,n+1,WAIT);
 }
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

 _DOUBLE_ dx = 1.0/n;
 double rp= kk*(b+1)*(b+1)/4;
 double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
 double dtr=1/(epsilon+((M1/M2)*rp));
 double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
 _DOUBLE_ dt = (_DOUBLE_) ddt;
 _DOUBLE_ alpha = d*dt/(dx*dx);

 // End Initization of various simulation variables

 // Report various information
 // Do not remove this call, it is needed for grading
if(rank==0)
  ReportStart(logfile, dt, T, m, n, tx, ty, noComm);

 // Start the timer
 double t0 = -MPI_Wtime();
 int niter = solve(logfile, &E, &E_prev, R, my_m, my_n, T, alpha, dt, do_stats, plot_freq,STATS_FREQ, &tmp_entire, rank,size, n,tx,ty,noComm);

 t0 += MPI_Wtime();

 // Report various information
 // Do not remove this call, it is needed for grading
 if(rank==0){
   if(noComm == 0)
     ReportEnd(logfile,T,niter,tmp_entire,m,n,t0,tx,ty);
   else
     ReportEnd(logfile,T,niter,E_prev,m,n,t0,tx,ty);

 }
 if (plot_freq){
    printf("\n\nEnter any input to close the program and the plot...");
    int resp;
    scanf("%d",&resp);
  }

 free (E);
 free (E_prev);
 free (R);

 MPI_Finalize();
}
