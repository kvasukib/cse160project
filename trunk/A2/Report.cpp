// 
// Performs various reporting functions
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "types.h"
using namespace std;

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem

 _DOUBLE_ stats(_DOUBLE_ **E, int m, int n, _DOUBLE_ *_mx);


// Report statistics periodically
void repNorms(ofstream& logfile,_DOUBLE_ **E,_DOUBLE_ t, _DOUBLE_ dt, int m,int n, int niter, int stats_freq){
     int k = (int)(t/stats_freq);
     if ((t-k*100)<dt) {
          _DOUBLE_ mx;
          _DOUBLE_ l2norm = stats(E,m,n,&mx);
          cout <<      setw(6);
          cout.setf(ios::fixed);
          cout << "iteration= " << niter << ", t= ";
          cout.unsetf(ios::fixed);
          cout.setf(ios::scientific);
          cout.precision(6);
          cout <<  t << endl;
          cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
          logfile <<          setw(6);
          logfile.setf(ios::fixed);
          logfile << "iteration= " << niter << ", t= ";
          logfile.unsetf(ios::fixed);
          logfile.setf(ios::scientific);
          logfile.precision(6);
          logfile <<  t << endl;
          logfile << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
     }
}
void printTOD(ofstream& logfile,string mesg)
{
        time_t tim = time(NULL);
        string s = ctime(&tim);
        if (mesg.length() ==  0) {
            cout << "Time of day: " << s.substr(0,s.length()-1) << endl;
            logfile << "Time of day: " << s.substr(0,s.length()-1) << endl;
        }
        else {
            cout << "[" << mesg << "] " ;
            cout << s.substr(0,s.length()-1) << endl;
            logfile << "[" << mesg << "] " ;
            logfile << s.substr(0,s.length()-1) << endl;
        }
        cout << endl;
        logfile << endl;
}


// Computes the gigaflops rate

double gflops(int n, int niter, double time){
    int n2 = n*n;
    int64_t updates = (int64_t) n2 * (int64_t) niter;
    int64_t flops = 28 * updates;
    double flop_rate = (double) flops / time;
    return ( flop_rate/1.0e9);
}


void ReportEnd(ofstream& logfile, _DOUBLE_ T,int niter, _DOUBLE_ **E_prev, int m,int n, double t0, int tx, int ty){
    _DOUBLE_ mx, l2norm;
    printTOD(logfile,"Simulation completes");
    l2norm = stats(E_prev,m,n,&mx);

    double gf = gflops(n+1, niter, t0);
    cout << "End at time: " << T;
    cout <<          setw(6);
    cout.setf(ios::fixed);
    cout << ", iteration " << niter << endl;
    cout.unsetf(ios::fixed);
    cout.setf(ios::scientific);
    cout.precision(5);
    cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
    cout.unsetf(ios::scientific);
    cout.unsetf(ios::fixed);
    cout.precision(6);
    cout << "Running Time: " << t0 << " sec." << endl;
    cout.precision(3);
    cout << "GFlop rate: " << gf << endl;;

    logfile << "End at time: " << T;
    logfile <<          setw(6);
    logfile.setf(ios::fixed);
    logfile << ", iteration " << niter << endl;
    logfile.unsetf(ios::fixed);
    logfile.setf(ios::scientific);
    logfile.precision(5);
    logfile << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
    logfile.unsetf(ios::scientific);
    logfile.unsetf(ios::fixed);
    logfile.precision(6);
    logfile << "Running Time: " << t0 << " sec." << endl;
    logfile.precision(3);
    logfile << "GFlop rate: " << gf << endl;
    logfile.precision(5);

    logfile << "   M x N, tx x ty, SimT, ";
    logfile <<  "#iter, T_p, Gflops, Linf, L2" << endl;
    logfile << "># " << m << " " << n << " ";
    logfile.precision(3);
    logfile << tx << " " << ty << " ";
    logfile.precision(6);
    logfile << T << " ";
    logfile << niter << " ";
    logfile.precision(4);
    logfile << t0 << " "  << gf << " ";

    logfile.unsetf(ios::fixed);
    logfile.setf(ios::scientific);
    logfile.precision(5);
    logfile << mx << " " << l2norm << endl;


}

void ReportStart(ofstream& logfile,_DOUBLE_ dt, _DOUBLE_ T, int m, int n, int tx, int ty){
  
#ifdef FLOAT
    cout << "Using Single precision arithmetic" << endl;
    logfile << "Using Single precision arithmetic" << endl;
#else
    cout << "Using Double precision arithmetic" << endl;
    logfile << "Using Double precision arithmetic" << endl;
#endif

    cout << "dt= " << dt << ", T = " << T << endl;
    cout << "The code will run for approximately " << (int)ceil(T/dt) << " timesteps" << endl;
    cout << "m x n = " << m << " x " << n << endl;
    cout << "thread geometry: " << tx << " x " << ty << endl;
    cout << endl;
    logfile << "dt= " << dt << ", T = " << T << endl;
    logfile << "The code will run for approximately " << (int)ceil(T/dt) << " timesteps" << endl;
    logfile << "m x n = " << m << " x " << n << endl;
    logfile << "thread geometry: " << tx << " x " << ty << endl;
    logfile << endl;
}
