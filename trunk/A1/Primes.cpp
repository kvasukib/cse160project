#include <stdlib.h>
#include <fstream>
#include <iostream>
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

// Globals
int64_t *candidates;
int *primes;
int n;



//
// --- Start of main module
//
int main(int argc, char *argv[]) {
    initInput(argc,argv);

    double t0 = -getTime();
    int i;
    int nPrimes = 0;
    for (i=0;i<n; i++) {
        if (isPrime (candidates[i])){
            nPrimes++;
            // If candidates[i] is prime, then set prime[i] to true
            primes[i] = TRUE;
        }
    }
    t0 += getTime();

    // The log file
    // Do not change the file name, as the autograder
    // assumes this name
    ofstream logfile("Log.txt",ios::out);
    ReportPrimes(logfile, nPrimes);
    ReportTimings(logfile,n-1,t0);

    cout << endl;
    delete [] candidates;
}
