/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "types.h"
#include "mpi.h"
using namespace std;

void repNorms(ofstream& logfile, _DOUBLE_ **E, _DOUBLE_ t, _DOUBLE_ dt, int m,int n, int niter, int stats_freq);


// Reports statistics about the computation: the L2 Norm and the Infinity NOrm
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
// The Linf norm is simply the maximum (absolute) value over
// all points in the array

 _DOUBLE_ stats(_DOUBLE_ **E, int m, int n, _DOUBLE_ *_mx){
     _DOUBLE_ mx = -1;
     _DOUBLE_ l2norm = 0;
     int i, j;
     for (j=1; j<=m+1; j++)
       for (i=1; i<=n+1; i++) {
	   l2norm += E[j][i]*E[j][i];
	   _DOUBLE_ fe = fabs(E[j][i]);
	   if (fe > mx)
	       mx = fe;
      }

     l2norm /= (_DOUBLE_) ((m+1)*(n+1));
     l2norm = sqrt(l2norm);

     *_mx = mx;
     return l2norm;
 }

int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, _DOUBLE_ T, _DOUBLE_ alpha, _DOUBLE_ dt, int do_stats, int plot_freq, int stats_freq, _DOUBLE_ *** _tmp_entire, int rank, int size, int full_n, int tx, int ty){


 // Simulated time is different from the integer timestep number
 _DOUBLE_ t = 0.0;
 // Integer timestep number
 int niter=0;

 _DOUBLE_ **E = *_E, **E_prev = *_E_prev;
 _DOUBLE_ ** tmp_entire;
   tmp_entire = *_tmp_entire;

 int requests;
 if(rank ==0 || rank == size -1)
    requests = 1;
 else
    requests = 2;
 MPI_Request* request_arr = (MPI_Request*)malloc(sizeof(MPI_Request)*requests);
 MPI_Status* status_arr = (MPI_Status*)malloc(sizeof(MPI_Status)*requests);
 MPI_Request send_request;
 //MPI_Request send_request, recv_request;
 //MPI_Status status;

int itertest=0;
 // We continue to sweep over the mesh until the simulation has reached
 // the desired simulation Time
 // This is different from the number of iterations
  while (t<T) {
/*  
#ifdef DEBUG
   printMat(E_prev,m,n);
   repNorms(logfile,E_prev,t,dt,m,n,niter, stats_freq);
   if (plot_freq)
    splot(E_prev,t,niter,m+1,n+1,WAIT);
#endif
*/
   t += dt;
   niter++;
   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    */
    

   int i,j;


/*
if(rank==1){
for(j=1; j<=m+1; j++)
{
    for(i=1; i<=n+1; i++)
        cerr <<  E_prev[j][i] << ' ';
    cerr << '\n';

}
cerr << "---------------------------------------------------------------\n";
}*/
if(tx == 1  || rank==0){
     for (j=1; j<=m+1; j++) 
       E_prev[j][0] = E_prev[j][2];
}else{//TODO

 }

if(tx==1 || rank == size-1){

      for (j=1; j<=m+1; j++) 
	E_prev[j][n+2] = E_prev[j][n];

}else{//TODO

 }
  //send ghost cells
  if(ty == 1 || rank==0)
  {
   for (i=1; i<=n+1; i++) 
      E_prev[0][i] = E_prev[2][i];
}else{
    MPI_Isend(E_prev[1],n+3,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, &send_request);
    MPI_Irecv(E_prev[0],n+3,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,request_arr);
  }
if(ty == 1 || rank== size-1)
  {
   for (i=1; i<=n+1; i++) 
      E_prev[m+2][i] = E_prev[m][i];

  }
  else
  {

    MPI_Isend(E_prev[m+1],n+3,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD, &send_request);
    if(rank==0)
      MPI_Irecv(E_prev[m+2],n+3,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,request_arr);
    else
      MPI_Irecv(E_prev[m+2],n+3,MPI_DOUBLE,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,request_arr+1);
  }

MPI_Waitall(requests, request_arr, status_arr); 
   // Solve for the excitation, a PDE
   for (j=1; j<=m+1; j++){
     for (i=1; i<=n+1; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
    }

   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   for (j=1; j<=m+1; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (i=1; i<=n+1; i++, EE++,RR++)
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }

   for (j=1; j<=m+1; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (i=1; i<=n+1; i++, EE++, RR++)
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
   }

/*   if(rank==0)
     for(i=1; i<= m+1; i++)
         for(j=1; j <= n+1; j++)
             E_prev[i][j] = 5;*/
/*   if(rank==1)
     for(i=1; i<= m+1; i++)
         for(j=1; j <= n+1; j++)
             E_prev[i][j] = 6;
*//*
   if(rank==0)
     for(i=0; i <=m+2; i++){
       for(j=0; j<=n+2; j++)
             cerr << E_prev[i][j];
       cerr << '\n';
     }
*/
/*
   if(rank==0){
     for(i = 0; i <=full_n+2; i++)
     {
        for(j=0; j <=full_n+2; j++){
            //printf("Reading i=%d j=%d\n inter=%f", i, j, t);
            cerr << tmp_entire[i][j];}
        cerr << '\n';
     }

   }*/
   //printf("POINTER VALUES %d %d\n",&E_prev[0][0],&tmp_entire[0][0]);
   //printf("VALUE IS %d   %d   %d\n",m+1,n+3, (m+1)*(n+3));
   int rc = MPI_Gather(&E[1][0],(m+1)*(n+3), MPI_DOUBLE, &tmp_entire[1][0],(m+1)*(n+3),MPI_DOUBLE,0,MPI_COMM_WORLD);
 /*
   if(rank==1)
     for(i=0; i <=m+2; i++){
       for(j=0; j<=n+2; j++)
             cerr << E_prev[i][j];
       cerr << '\n';
     }
  */
   if(rank==0){
     if (do_stats)
       repNorms(logfile, tmp_entire,t,dt,full_n,full_n,niter,stats_freq);
   
     /*for(i = 0; i <=full_n+2; i++)
     {
        for(j=0; j <=full_n+2; j++){
            //printf("Reading i=%d j=%d\n", i, j);
            printf("%1.1f ",tmp_entire[i][j]);}
        printf("\n");

     }*/
     if (plot_freq){
          int k = (int)(t/plot_freq);
          if ((t-k*plot_freq)<dt){
              splot(tmp_entire,t,niter,full_n+1,full_n+1,WAIT);
          }
      }
   }
   //printf("THREAD %d GOT HERE\n", rank);
   // Swap current and previous
   _DOUBLE_ **tmp = E; E = E_prev; E_prev = tmp;
   //char c;
   //cin >> c;
 }

  // Store them into the pointers passed in
  *_E = E;
  *_E_prev = E_prev;
  
  return niter;
}
