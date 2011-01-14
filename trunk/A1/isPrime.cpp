#include <stdlib.h>
#include <iostream>
using namespace std;
//
// Tests for primality
// Based on Java code written by Alan Kaminsky
// Ported to C++ by Scott B. Baden, UCSD, Dec. 28, 2010
//

#define FALSE 0
#define TRUE 1

//
// Returns true if x is prime, else false
int isPrime(int64_t x) {
    if (x == 2) return TRUE;
    if (x % 2 == 0) return FALSE;
    int64_t p = 3;
    int64_t psqr = p*p;
    while (psqr <= x) {
        if (x % p == 0) return FALSE;
        p += 2;
        psqr = p*p;
    }
    return TRUE;
}
