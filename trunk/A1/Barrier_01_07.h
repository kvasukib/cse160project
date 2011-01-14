#ifndef _BARRIER_H
#define _BARRIER_H

/****************************************************************
* Barrier.h							*
*								*
* Barrier Synchronization for POSIX threads			*
*								*
* Author:  Pietro Cicotti 		*
*								*
*****************************************************************/

#include <cassert>
#include "Semaphore.h"

class Barrier
{
    int count;
    int capacity;
    Semaphore arrivals;
    Semaphore departures;

public:

    Barrier(int threads =2) : arrivals(1), departures(0) {
        count = 0;
        capacity = threads;
    }

    void synchronize() {
        arrivals.P();
        ++count; // atomically count the witing threads

        if(count < capacity)
            arrivals.P();
        else // last  processor enables all to go	
            departures.V();

        departures.P();
        --count;  // atomically decrement
        if(count > 0)
            departures.V();
        else
            arrivals.V();
    }

};

#endif
