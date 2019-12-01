#include <iostream>
#include <sstream>
#include <vector>
#include <ctime>
#include "mylib.h"
using namespace std;

int main(int argc, char **argv)
{
    int const NLOOPS = 50;        // chose a value such that the benchmark runs at least 10 sec.
    long long int N = 40000000;
    vector<double> a(N);

	// vector initialiation
	for (int i=1;i<=N;i++)
	{
		a[i] = (1.0/i)*(1.0/i);
	}

    double tstart, t1;  // timer

    // Do calculation
    tstart = clock();       // start timer
 
    for (int i=0;i<NLOOPS;i++)
    {
		Kahan_scalar(a,N);
	}
    t1 = clock() - tstart;
  
//##########################################################################
// Timings  and Performance  
    t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
    cout.precision(6);
    cout << "Total time in sec: " << t1 << endl;
    t1 /= NLOOPS; // divide by number of function calls
    cout.precision(6);
    cout << "Timing in sec. : " << t1 << endl;
    cout << "Memory allocated (in GB): " << 2.0*N*sizeof(double)/1024/1024/1024 << endl;
    cout << "GFLOPS: " << 4.0*N/t1/1024/1024/1024 << endl;
    cout << "Gib/sec: "<< 12.0*N*sizeof(double)/t1/1024/1024/1024 << endl;
//########################################################################## 
    cout << "Normal summation result is:" << scalar(a,N)<< endl;
    cout << "Kahan summation result is:"  << Kahan_scalar(a,N)<< endl;
    cout << "The difference: KahanSum - NormalSum = " << Kahan_scalar(a,N) - scalar(a,N)<<endl;
    cout << "The difference: pi^2/6 - NormalSum = " << pow(M_PI,2.0)/6.0 - scalar(a,N)<< endl;
    cout << "The difference: pi^2/6 - KahanSum = " << pow(M_PI,2.0)/6.0 - Kahan_scalar(a,N)<<endl;

    return 0;
}
