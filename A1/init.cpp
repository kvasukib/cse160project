// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <stdlib.h>
#include <assert.h>
#include <iostream>

using namespace std;
// Input the data

#define FALSE 0
#define TRUE 1

// Shared globals, allocated before threads are spawned
extern int64_t *candidates;
extern int *primes;
extern int n;
void initInput(int argc, char *argv[]){
    n = argc; n--;
    cout << "Testing " << n << " candidate primes" << endl;

// Shared global, allocated before threads are spawned
    candidates = new int64_t[n];
    int i;
    for (i=0;i<n; i++) {
        candidates[i] = atoll(*++argv);
        cout << "   Candidate prime # " << i << " is " << candidates[i] << endl;
    }

    primes = new int[n];
    for (i=0;i<n; i++)
        primes[i] = FALSE;
}
