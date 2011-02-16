// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
// A parallel timer
// Times the execution of the whole program, using  a barrier
// to sycnhronize the threads
//
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <pthread.h>
#ifndef PTHREAD_BARRIER
// Use our own barrier if the local pthreads implementation doesn't define one
#include "Barrier.h"
#endif

using namespace std;

const double kMicro = 1.0e-6;
#ifdef  PTHREAD_BARRIER
pthread_barrier_t barrTimer;
#else
Barrier* barrTimer;
#endif


// Prior to using the timer, you must invoke initTimer() exactly
// one time prior to spawning threads
// We should make a timer package and make this cleaner
// to avoid the need for initialization

void initTimer(int NT)
{
#ifdef PTHREAD_BARRIER
    if(pthread_barrier_init(&barrTimer, NULL, NT)) {
        cerr << "Could not create a barrier\n";
        exit(-1); // return -1;
    }
#else
    barrTimer = new Barrier(NT);
#endif
}

// Before the program terminates invoke FinalizeTimer() exactly
// one time after joining all the spawned threads
void FinalizeTimer()
{
#ifdef PTHREAD_BARRIER
    pthread_barrier_destroy(&barrTimer);
#else
    delete barrTimer;
#endif
}

double getTime(int DUMMY)
{
    struct timeval TV;

    // We should really test the return code, but
    // that the standard doesn't specify precisely
    // which threads is returned which value!
#ifdef PTHREAD_BARRIER
    int rc = pthread_barrier_wait(&barrTimer);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
        cerr << "Could not wait on barrier\n";
        exit(-1); // return -1;
    }
#else
    barrTimer->synchronize();
#endif

    const int RC = gettimeofday(&TV, NULL);
    if(RC == -1)
    {
        cerr << "ERROR: Bad call to gettimeofday()\n";
        return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()
