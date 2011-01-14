#ifndef _SEMAPHORE_WRAPPER_H
#define _SEMAPHORE_WRAPPER_H

/****************************************************************
* Semaphore.h							*
*								*
* Semphore class, a wrapper to <semaphore.h> 			*
*								*
* Author:  Pietro Cicotti 					*
*  Modified by Scott Baden to work with semaphore.h		*
*								*
*****************************************************************/

#include <cassert>
#include <semaphore.h>

class Semaphore
{
   sem_t sem;

public:

   Semaphore(int c =0) {
       // Initialize the count to c
       // 2nd parameter should be 0
       sem_init(&sem,0,c);
   }

#if 0
   ~Semaphore() {
   }

#endif
   void P() { 
	sem_post(&sem);
   }

   void V() { 
	sem_wait(&sem);
   }

};

#endif
