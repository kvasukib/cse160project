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
#include <pthread.h>
#include "time.h"
#include "apf.h"
#include "types.h"
using namespace std;

extern pthread_barrier_t barr;
//global vars
extern int q, first_q, rest_q; //used to determine partitioning
extern _DOUBLE_ **E, **R, **E_prev;
 extern _DOUBLE_ T;
 extern int m;
extern  int n;
extern  int do_stats;
extern  int plot_freq;
extern  int tx, ty;

extern  _DOUBLE_ dt;
extern  _DOUBLE_ alpha;

// This parameter controls the frequncy (in timesteps)
// that summary statistics are reported
// Reduce the value of FREQ to increase the frequency,
// increase the value to raise the frequency
const int STATS_FREQ = 100;

extern ofstream logfile;

extern int NT;


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

 // Simulated time is different from the integer timestep number
 _DOUBLE_ t = 0.0;
 // Integer timestep number
extern int niter;

//int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, _DOUBLE_ T, _DOUBLE_ alpha, _DOUBLE_ dt, int do_stats, int plot_freq, int stats_freq){
void * solve_thr (void * arg){
 

 //_DOUBLE_ **E = *_E, **E_prev = *_E_prev;
  int64_t id = reinterpret_cast<int64_t>(arg);
  int tid = id;
 
  int start;
  int end;
  if(tid < q)
  {
    start = 1 + id * first_q;
    end = start + first_q;
  }
  else
  {
    start = 1+ (q*first_q) + (id-q)*rest_q;
    end = start+ rest_q;
  }
 // We continue to sweep over the mesh until the simulation has reached
 // the desired simulation Time
 // This is different from the number of iterations
  while (t<T) {
  //int asd = 0;
  //cerr << " ";
  if(tid == 0){  
   #ifdef DEBUG
   printMat(E_prev,m,n);
   repNorms(logfile,E_prev,t,dt,m,n,niter, STATS_FREQ);
   if (plot_freq)
    splot(E_prev,t,niter,m+1,n+1,WAIT);
    #endif


    niter++;

   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    */
    
     int i,j;
     for (j=1; j<=m+1; j++) 
       E_prev[j][0] = E_prev[j][2];

      for (j=1; j<=m+1; j++) 
	E_prev[j][n+2] = E_prev[j][n];
 
   for (i=1; i<=n+1; i++) 
      E_prev[0][i] = E_prev[2][i];

   for (i=1; i<=n+1; i++) 
      E_prev[m+2][i] = E_prev[m][i];
   }//end first serial section

   int rc = pthread_barrier_wait(&barr);
   if (rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
   {
      cerr << "could not wait on barrier\n";
      exit(-1);
   }

   
   // Solve for the excitation, a PDE
   for (int j=start; j<end; j++){
     for (int i=1; i<=n+1; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
    }


   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   for ( int j=start; j<end; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (int i=1; i<=n+1; i++, EE++,RR++)
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }

   for (int j=start; j<end; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (int i=1; i<=n+1; i++, EE++, RR++)
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
   }
    if(tid == 0)
      t += dt;
   //barrier to finish all computations
   rc = pthread_barrier_wait(&barr);
   if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
   {
      cerr << "could not wait on barrier\n";
      exit(-1);
   }

   //cerr << tid;
   if(tid == 0)
   {
     if (do_stats)
       repNorms(logfile, E,t,dt,m,n,niter,STATS_FREQ);

     if (plot_freq){
          int k = (int)(t/plot_freq);
          if ((t-k*plot_freq)<dt){
            splot(E,t,niter,m+1,n+1,WAIT);
          }
      }

   // Swap current and previous
     _DOUBLE_ **tmp = E; E = E_prev; E_prev = tmp;

   }//end serial section


 }
 //pthread_barrier_destroy(&barr);
 //cerr << "thread "<< tid <<" exiting\n";
  // Store them into the pointers passed in
  //*_E = E;
  //*_E_prev = E_prev;
  pthread_exit(NULL);
  return 0;
}
