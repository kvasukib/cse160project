// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
const double kMicro = 1.0e-6;

double getTime()
{
    struct timeval TV;

    const int RC = gettimeofday(&TV, NULL);
    if(RC == -1)
    {
        cerr << "ERROR: Bad call to gettimeofday()\n";
        return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()
