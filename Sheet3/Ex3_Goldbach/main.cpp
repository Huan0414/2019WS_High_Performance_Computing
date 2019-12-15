#include <algorithm> 
#include <vector>
#include "mylib.h"
#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <sstream>


using namespace std;


int main()
{
    int n;
    cout << "Please enter a even number larger than 3."<<"\n";
    cin  >> n;
    cout << endl;

//########################################################################
    //int nthreads;                                  // OpenMP
    //omp_set_num_threads(64);
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

//#################################################################### 
//all numbers between 4-n decompositions conditions
    double t1, tstart; //time measurement
    
	int N = (n-4)/2 + 1;        //Only consider even numbers from 4-n
    vector<int> pair_cont(N,0); // all counts are initialized as zero
    auto primes = get_primes(n); // get all the primes until n
    
    tstart  = omp_get_wtime();
    count_goldbach(n, pair_cont, primes); 
    
    t1 = omp_get_wtime() - tstart;
    //t1 /= CLOCKS_PER_SEC;
       

//#################################################################### 
//Check the number with maximum pairs
    cout << "The even number has most prime pairs is" <<" "<< 2*(2 + distance(pair_cont.begin(), max_element(pair_cont.begin(),pair_cont.end())))<< endl;
	cout << "It has " << *max_element(pair_cont.begin(),pair_cont.end()) << " prime pairs" << endl;
//####################################################################    
// Performance evaluation
    cout << endl;
    cout.precision(4);
    cout << "Total time	:" << t1 << endl;
    //cout << "Time per loop	:" << t2 << endl;
	
    return 0;
}
