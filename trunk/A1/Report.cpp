// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iostream>

using namespace std;

#define FALSE 0
#define TRUE 1


// Shared globals, allocated before threads are spawned
extern int64_t *candidates;
extern int *primes;
extern int n;


void ReportPrimes(ofstream& logfile, int nPrimes){
    int i;
    int np = 0;
    // Check to see if all numbers are prime.
    // If so, then we output that message,
    // but without cluttering the output with a message for every 
    // prime we tested.
    // If not all numbers were prime, then we rescan and do a full printout
    // of those numbers that were prime
    // Most of this output also goes to the logfile

    for (i=0;i<n; i++) {
        if (primes[i]){
            np++;
        }
    }
    if (np != nPrimes){
        cerr << "Sanity check failure" << endl;
    }
    if (nPrimes == n){
        cout << endl <<  "*** All " << n << " numbers were prime" << endl;
        logfile << endl <<  "*** All " << n << " numbers were prime" << endl;
        logfile << "#: " << nPrimes << endl;
    }
    else{
        cout << endl;
        logfile << endl;
        cout << "A total of " << nPrimes << " primes were identified" << endl;
        logfile << "#: " << nPrimes << endl;
        for (i=0;i<n; i++) {
            if (primes[i]){
                cout << candidates[i] << " is a prime number" << endl;
                logfile << candidates[i] << " is a prime number" << endl;
            }
        }
    }
}

void ReportTimings(ofstream& logfile, int n, double timings){
    cout << "Testing the " << n << " candidates " << " took ";
    cout << "  " << timings << " seconds" << endl;
    logfile << "T: " << timings << endl;
    }
