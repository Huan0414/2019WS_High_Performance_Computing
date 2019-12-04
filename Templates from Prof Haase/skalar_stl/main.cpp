#include <iostream>
#include <sstream>

#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <ctime>
#include "mylib.h"
using namespace std;

int main(int argc, char **argv)
{
    int const NLOOPS = 50;        // chose a value such that the benchmark runs at least 10 sec.
    unsigned int N = 50000001;
//##########################################################################
//   Read Paramater from command line  (C++ style)
    cout << "Checking command line parameters for: -n <number> " << endl;
    for (int i = 1; i < argc; i++)
    {
        cout << " arg[" << i << "] = " << argv[i] << endl;
        if (std::strncmp(argv[i], "-n", 2) == 0 && i + 1 < argc) // found "-n" followed by another parameter
        {
            N = static_cast<unsigned int>(atoi(argv[i + 1]));
        }
        else
        {
            cout << "Corect call: " << argv[0] << " -n  <number>\n";
        }
    }

    cout << "\nN = " << N << endl;

//##########################################################################
//  Memory allocation
    cout << "Memory allocation\n";

    vector<double> x(N), y(N);

    cout.precision(2);
    cout << 2.0 * N *sizeof(x[0]) / 1024 / 1024 / 1024 << " GByte Memory allocated\n";
    cout.precision(6);

//##########################################################################
//  Data initialization
//  Special:  x_i = i+1;  y_i = 1/x_i  ==> <x,y> == N
    for (unsigned int i = 0; i < N; ++i)
    {
        x[i] = i + 1;
        y[i] = 1.0 / x[i];
    }

//##########################################################################
    cout << "\nStart Benchmarking\n";

    double tstart, t1;  // timer
    double sk(0.0);

// Do calculation
    tstart = clock();       // start timer

    for (int i = 0; i < NLOOPS; ++i)
    {
        sk = scalar(x, y);
//      sk = norm(x);
    }

    t1 = clock() - tstart;  // stop timer
    t1 /= CLOCKS_PER_SEC;   // now, t1 in seconds
    t1 /= NLOOPS;           // divide by number of function calls

//##########################################################################
// Check the correct result
    cout << "\n <x,y> = " << sk << endl;
    if (static_cast<unsigned int>(sk) != N)
    {
        cout << "  !!   W R O N G  result   !!\n";
    }
    cout << endl;

//##########################################################################
// Timings  and Performance
    cout << endl;
    cout.precision(2);
    cout << "Timing in sec. : " << t1 << endl;
    cout << "GFLOPS         : " << 2.0 * N / t1 / 1024 / 1024 / 1024 << endl;
    cout << "GiByte/s        : " << 2.0 * N / t1 / 1024 / 1024 / 1024 * sizeof(x[0]) << endl;

//##########################################################################
    return 0;
}  // memory for x and y will be deallocated their destructors
