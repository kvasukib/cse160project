#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <pthread.h>

#include <math.h>

using namespace std;
// Tests numbers for primality
// Usage: a.out n1 [n2...]
// Where the n_i are numbers to be tested for primality
// Based on Java code written by Alan Kaminsky
// Ported to C++ by Scott B. Baden, UCSD, Dec. 28, 2010
//

#define FALSE 0
#define TRUE 1

int isPrime(int64_t x);
double getTime();

void initInput(int argc, char *argv[]);
void ReportPrimes(ofstream& logfile, int nPrimes);
void ReportTimings(ofstream& logfile, int n, double timings);

void initTimer(int NT);
void FinalizeTimer();
void *prime_thr( void *arg );

// Globals
int64_t *candidates;
int *primes;
int n;
int NT;

//our globals
int elements_per_thread;
double t0;

//
// --- Start of main module
//
int main(int argc, char *argv[]) {
    initInput(argc,argv);

    int i;
    int nPrimes = 0;
    pthread_t * th_arr = new pthread_t [NT];
    elements_per_thread = (int) ceil(n/NT);

    initTimer(NT);
    //double t0 = -getTime();
    for (i=0;i<NT; i++) {
        int64_t ind = i;
        pthread_create(&th_arr[i], NULL, prime_thr, reinterpret_cast<void *>(ind));

    }
    
    // Join the threads
    for (int t=0; t<NT; t++)
      pthread_join(th_arr[t],NULL);

    FinalizeTimer();
    //t0 += getTime();
    for (int i=0; i<n; i++){
      if (primes[i] == TRUE){
        nPrimes ++;
      }
    }


    // The log file
    // Do not change the file name, as the autograder
    // assumes this name
    ofstream logfile("Log.txt",ios::out);
    ReportPrimes(logfile, nPrimes);
    ReportTimings(logfile,n-1,t0);

    cout << endl;
    delete [] candidates;

    pthread_exit(NULL);
}
