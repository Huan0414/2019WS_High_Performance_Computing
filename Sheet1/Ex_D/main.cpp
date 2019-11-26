#include <iostream>
#include <sstream>
#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <ctime>
#include "mylib.h"
using namespace std;

double scalar(long long int n);
double Kahan_scalar(long long int n);

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
    cout << "\nStart Benchmarking\n";

    double tstart, t1;  // timer

// Do calculation
    tstart = clock();       // start timer
    int n;
    bool infinity;

    cout << "Please Enter 0 or 1.\n"
         << "Enter 0 for normal case;\n"
         << "Enter 1 for infinite n."<<endl;
    cin  >> infinity;


    if (!infinity)
    {
    cout << "Please enter n."<<endl;
    cin  >> n;
    cout << "Normal summation result is:" << scalar(n)<< endl;
    cout << "Kahan summation result is:"  << Kahan_scalar(n) << endl;
    cout << "The difference: KahanSum - NormalSum = " << Kahan_scalar(n) - scalar(n)<<endl;
    }
    else
    {
    n = 9999999;
    cout << "Normal summation result is:" << scalar(n)<< endl;
    cout << "Kahan summation result is:"  << Kahan_scalar(n) << endl;
    cout << "The difference: KahanSum - NormalSum = " << Kahan_scalar(n) - scalar(n)<<endl;
    cout << "The difference: pi^2/6 - NormalSum = " << pow(M_PI,2.0)/6.0 - scalar(n)<< endl;
    cout << "The difference: pi^2/6 - KahanSum = " << pow(M_PI,2.0)/6.0 - Kahan_scalar(n)<<endl;
    }

    t1 = clock() - tstart;  // stop timer
    t1 /= CLOCKS_PER_SEC;   // now, t1 in seconds
    t1 /= NLOOPS;           // divide by number of function calls


//##########################################################################
// Timings  and Performance
    cout << endl;
    cout
    .precision(2);
    cout << "Timing in sec. : " << t1 << endl;
    cout << "GFLOPS         : " << 2.0 * N / t1 / 1024 / 1024 / 1024 << endl;
    cout << "GiByte/s        : " << 2.0 * N / t1 / 1024 / 1024 / 1024 * sizeof(x[0]) << endl;


//##########################################################################
    return 0;
}
