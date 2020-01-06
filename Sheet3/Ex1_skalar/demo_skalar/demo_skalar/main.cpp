#include "check_env.h"
#include "mylib.h"
#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <sstream>
#include <string>
using namespace std;

int main(int argc, char const *argv[])
{
	omp_set_num_threads(10);
	cout << "Number of available processors: " << omp_get_num_procs() << endl;
	cout << "The return of command 'omp_in_parallel': " << omp_in_parallel << endl;
    int const NLOOPS = 200;        // chose a value such that the benchmark runs at least 10 sec.
    unsigned int N = 40000000;
//##########################################################################
//   Read Parameter from command line  (C++ style)
    cout << "Checking command line parameters for: -n <number> " << endl;
    for (int i = 1; i < argc; i++)
    {
        cout << " arg[" << i << "] = " << argv[i] << endl;
        string ss(argv[i]);
        if ("-n"==ss && i + 1 < argc) // found "-n" followed by another parameter
        {
            N = static_cast<unsigned int>(atoi(argv[i + 1]));
        }
        else
        {
            cout << "Corect call: " << argv[0] << " -n  <number>\n";
        }
    }

    cout << "\nN = " << N << endl;
    
    check_env(argc, argv);
////########################################################################
    //int nthreads;                                  // OpenMP
    //omp_set_num_threads(10);
    //#pragma omp parallel default(none) shared(cout,nthreads)
    //{
        //int const th_id  = omp_get_thread_num();   // OpenMP
        //int const nthrds = omp_get_num_threads();  // OpenMP
        
        //stringstream ss;
        //ss << "C++: Hello World from thread " << th_id << " / " << nthrds << endl;
        //#pragma omp critical
        //{
            //cout << ss.str();                      // output to a shared ressource
        //}
        //#pragma omp master
        //nthreads = nthrds;                         // transfer nn to to master thread
    //}
    //cout << "   " << nthreads << "   threads have been started." << endl;

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
//	scalarNofor(x,y);
//	cout << "scalarNoFor: " << sum << endl; 
    
// Do calculation
    double tstart = omp_get_wtime();                  // OpenMP

    double sk(0.0);
    for (int i = 0; i < NLOOPS; ++i)
    {
//        sk = scalar(x, y);
        sk = scalarNofor(x,y);
//      sk = norm(x);
    }

    double t1 = omp_get_wtime() - tstart;             // OpenMP
    double t2 = t1/NLOOPS;           // divide by number of function calls

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
    cout << "Total time in sec.	: " << t1 << endl;
    cout << "Time/loop in sec.	: " << t2 << endl;
    cout << "GFLOPS		        : " << 2.0 * N / t2 / 1024 / 1024 / 1024 << endl;
    cout << "GiByte/s	        : " << 2.0 * N / t2 / 1024 / 1024 / 1024 * sizeof(x[0]) << endl;

//#########################################################################

    cout << "\n  Try the reduction with an STL-vektor \n";
    
    auto vr = reduction_vec(100);
    cout << "done\n";
    //cout << vr << endl;
    
    auto vra1 = reduction_vec_append(3);
    cout << "Appending result is: " << vra1 << endl;
    
    auto vra2 = reduction_vec_append_manually(3);
    cout << "Manually ordered appending result is: " << vra2 << endl;
    
    //vector<int> vr1(3);
    //vr1 = {1,2,3};
    //vector<int> vr2 = {4,5,6};
    //cout << "Print vr1 : " << vr1 << endl;
    //cout << "Print vr2 : " << vr2 << endl;
    
    //vr1.insert(vr1.end(),vr2.begin(),vr2.end());
    
    //cout << "Apending vr1 & vr2 : " << vr1 << endl;

    return 0;
    
}  // memory for x and y will be deallocated their destructors
