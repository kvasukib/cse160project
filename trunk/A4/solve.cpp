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
_DOUBLE_ **alloc2D(int m,int n);

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

int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, _DOUBLE_ T, _DOUBLE_ alpha, _DOUBLE_ dt, int do_stats, int plot_freq, int stats_freq, _DOUBLE_ *** _tmp_entire, int rank, int size, int full_n, int tx, int ty, int noComm){


 // Simulated time is different from the integer timestep number
 _DOUBLE_ t = 0.0;
 // Integer timestep number
 int niter=0;

 _DOUBLE_ **E = *_E, **E_prev = *_E_prev;
 _DOUBLE_ ** tmp_entire;
   tmp_entire = *_tmp_entire;

 _DOUBLE_ * ghostsL, *ghostsR;
 if(ty == 1){//vertical partitioning ghosts
  if(rank!=0)
       ghostsL = (_DOUBLE_*)malloc(sizeof(_DOUBLE_)*(m+3));
  if(rank!=size-1)
       ghostsR = (_DOUBLE_*)malloc(sizeof(_DOUBLE_)*(m+3));

}
 int nrequests;
 if(rank ==0 || rank == size -1)
    nrequests = 1;
 else
    nrequests = 2;
 MPI_Request* request_arr = (MPI_Request*)malloc(sizeof(MPI_Request)*nrequests);
 MPI_Status* status_arr = (MPI_Status*)malloc(sizeof(MPI_Status)*nrequests);
 MPI_Request send_request;
 //MPI_Request send_request, recv_request;
 //MPI_Status status;

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
if(noComm || tx == 1  || rank==0){
     for (j=1; j<=m+1; j++) 
       E_prev[j][0] = E_prev[j][2];
}else{
    for (j = 0;j <= m+2;j++)
        ghostsL[j] = E_prev[j][1];
    MPI_Isend(ghostsL,m+3,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&send_request);
    MPI_Irecv(ghostsL,m+3,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,request_arr);
 }

if(noComm || tx==1 || rank == size-1){

      for (j=1; j<=m+1; j++) 
	E_prev[j][n+2] = E_prev[j][n];

}else{
    for (j = 0;j <= m+2;j++)
        ghostsR[j] = E_prev[j][n+1];
    MPI_Isend(ghostsR,m+3,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&send_request);
    if (rank == 0){
        MPI_Irecv(ghostsR,m+3,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,request_arr);
    }
    else{
        MPI_Irecv(ghostsR,m+3,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,request_arr+1);
    }

 }
  //send ghost cells
  if(noComm || ty == 1 || rank==0)
  {
   for (i=1; i<=n+1; i++) 
      E_prev[0][i] = E_prev[2][i];
}else{
    MPI_Isend(E_prev[1],n+3,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, &send_request);
    MPI_Irecv(E_prev[0],n+3,MPI_DOUBLE,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,request_arr);
  }
if(noComm || ty == 1 || rank== size-1)
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

   // Solve for the excitation, a PDE
   for (j=2; j<=m; j++){
     for (i=2; i<=n; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
    }

   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   for (j=2; j<=m; j++){
     double *RR = &R[j][2];
     double *EE = &E[j][2];
     for (i=2; i<=n; i++, EE++,RR++)
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }

   for (j=2; j<=m; j++){
     double *RR = &R[j][2];
     double *EE = &E[j][2];
     for (i=2; i<=n; i++, EE++, RR++)
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
   }

if(noComm == 0)
  MPI_Waitall(nrequests, request_arr, status_arr); 

if (ty == 1 && rank != 0)
    for (j = 0;j <= m+2;j++)
        E_prev[j][0] = ghostsL[j];

if (ty == 1 && rank != size - 1)
    for (j = 0;j <= m+2;j++)
        E_prev[j][n+2] = ghostsR[j];


   // Solve for the excitation, a PDE
   for (j=1; j<=m+1; j+=m){
     for (i=1; i<=n+1; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
      }
   }
   // Solve for the excitation, a PDE
   for (j=2; j<=m; j++){
     for (i=1; i<=n+1; i+=n) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
   }



   for (j=1; j<=m+1; j+=m){
     double *RR = &R[j][1];
     double *EE = &E[j][1];
     for (i=1; i<=n+1; i++, EE++,RR++)
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }

   for (j=2; j<=m;j++){
     double *RR = &R[j][1];
     double *EE = &E[j][1];
     for (i=1; i<=n+1;i+=n,EE+=n,RR+=n)
        EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }
   for (j=1; j<=m+1; j+=m){
     double *RR = &R[j][1];
     double *EE = &E[j][1];
     for (i=1; i<=n+1; i++, EE++, RR++)
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
   }

   for (j=2; j<=m;j++){
     double *RR = &R[j][1];
     double *EE = &E[j][1];
     for (i=1;i<=n+1;i+=n,EE+=n,RR+=n)
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
   //int rc = MPI_Gather(&E[1][0],(m+1)*(n+3), MPI_DOUBLE, &tmp_entire[1][0],(m+1)*(n+3),MPI_DOUBLE,0,MPI_COMM_WORLD);
 /*
   if(rank==1)
     for(i=0; i <=m+2; i++){
       for(j=0; j<=n+2; j++)
             cerr << E_prev[i][j];
       cerr << '\n';
     }
  */
if(do_stats || plot_freq){
   if(tx == 1 && noComm == 0)
     int rc = MPI_Gather(&E[1][0],(m+1)*(n+3), MPI_DOUBLE, &tmp_entire[1][0],(m+1)*(n+3),MPI_DOUBLE,0,MPI_COMM_WORLD);
   else if (ty==1 && noComm == 0){
     //vertical gather
   }
   if(rank==0){
     if (do_stats){
       if(noComm ==0)
         repNorms(logfile, tmp_entire,t,dt,full_n,full_n,niter,stats_freq);
       else
         repNorms(logfile, E,t,dt,full_n,full_n,niter,stats_freq);
     }
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
              if(noComm == 0)
                splot(tmp_entire,t,niter,full_n+1,full_n+1,WAIT);
              else
                splot(E,t,niter,full_n+1,full_n+1,WAIT);
      
          }
      }
   }
}
   //printf("THREAD %d GOT HERE\n", rank);
   // Swap current and previous
   _DOUBLE_ **tmp = E; E = E_prev; E_prev = tmp;
   //char c;
   //cin >> c;
 }
if(tx==1 && noComm == 0)
   int rc = MPI_Gather(&E_prev[1][0],(m+1)*(n+3), MPI_DOUBLE, &tmp_entire[1][0],(m+1)*(n+3),MPI_DOUBLE,0,MPI_COMM_WORLD);
else if (ty==1 && noComm == 0)
{   //vertical gather

       _DOUBLE_ *tmp_block = (_DOUBLE_*)malloc(sizeof(_DOUBLE_)*(m+3)*(n+1));
       _DOUBLE_ * tmp_vert = (_DOUBLE_*)malloc(sizeof(_DOUBLE_)*(m+3)*(n+1)*size);
    
   int k = 0;
   for(int j = 0; j <=m+2; j++)
     for(int i = 1; i <= n+1; i++)
       tmp_block[k++] = E_prev[j][i];


   int rc = MPI_Gather(tmp_block,(m+3)*(n+1), MPI_DOUBLE, tmp_vert,(m+3)*(n+1),MPI_DOUBLE,0,MPI_COMM_WORLD);


        if (rank == 0){
            int source;
            int start;
            int end;
            int iter = 0;
            for (int i = 0;i < size;i++){
                source = i;
                start = (n+1)*source +1;
                end = start + n;
                for (int j = 0;j <=m+2;j++)
                    for (int i = start;i <=end;i++){
                        tmp_entire[j][i] = tmp_vert[iter++];
                    }
            }
        }
}  
/*
   for(int i = 0; i <=full_n+2; i++)
     {
        for(int j=0; j <=full_n+2; j++){
            //printf("Reading i=%d j=%d\n", i, j);
            printf("%1.1f ",tmp_entire[i][j]);}
        printf("\n");

     }*/
  // Store them into the pointers passed in
  *_E = E;
  *_E_prev = E_prev;
  
  return niter;
}
