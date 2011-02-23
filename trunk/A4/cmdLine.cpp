// Process command line arguments
// 
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "types.h"
using namespace std;

void cmdLine(int argc, char *argv[], _DOUBLE_& T, int& n, int& tx, int& ty, int& do_stats, int& plot_freq, int& noComm){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"tx", required_argument, 0, 'x'},
        {"ty", required_argument, 0, 'y'},
        {"tfinal", required_argument, 0, 't'},
        {"stats", no_argument, 0, 's'},
        {"plot", required_argument, 0, 'p'},
        {"nocomm", no_argument, 0, 'k'},
 };
    // Process command line arguments
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:x:y:i:j:t:skp:",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                n = atoi(optarg);
                break;

	    // X processor geometry
            case 'x':
                tx = atoi(optarg);
                break;

	    // X processor geometry
            case 'y':
                ty = atoi(optarg);
                break;

	    // Length of simulation, in simulated time units
            case 't':
                T = atof(optarg);
                break;

	    // Print various statistics
            case 's':
                do_stats = 1;
                break;

	    // Shut off communication
            case 'k':
                noComm = 1;
                break;


	    // Plot the excitation variable
            case 'p':
                plot_freq = atoi(optarg);
                break;

	    // Error
            default:
                cout << "Usage: a.out [-n <domain size>] [-t <final time >]";
                cout << "\n\t    ";
                cout << " [-s print statistics] [-p <plot frequency>]\n\t";
                cout << "     [-tx <x processor geometry> [-ty <y processor geometry] [-k <no communication>]" << endl;
                exit(-1);
            }
    }
 }
}
