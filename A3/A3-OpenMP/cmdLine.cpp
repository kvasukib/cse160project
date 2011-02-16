//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
// Process command line arguments
// 
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "lu.h"
using namespace std;

// Process command line arguments
void cmdLine(int argc, char *argv[], int& n, int& nt, int& matCode) {
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"nt", required_argument, 0, 't'},
        {"l", no_argument, 0, 'l'},
        {"h", no_argument, 0, 'h'},

 };

 // Default values
 n = 16; nt  = 1; matCode = D_MAT;
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:t:lh",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                n = atoi(optarg);
                break;

	    // # threads
            case 't':
                nt = atoi(optarg);
                break;

	    // Use the 'Laplacian' matrix (5-point Laplacian in 2D)
            // Doesn't need partial pivoting
            case 'l':
                matCode = L_MAT;
                break;

	    // Use a Hadamard matrix 
            // Requires partial pivoting
            case 'h':
                matCode = H_MAT;
                break;

	    // Error
            default:
                cout << "Usage: lu [-n <domain size>] [-t <# threads>] [-l <use 5-point Laplacian>] [-h <use Hadamard Matrix>]";
                cout << endl;
                exit(-1);
            }
    }
 }
}
