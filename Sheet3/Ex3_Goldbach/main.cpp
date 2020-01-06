#include <algorithm> 
#include <vector>
#include "mylib.h"
#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <sstream>


using namespace std;

//g++ -O3 -fopenmp main.cpp mylib.h get_primes.cpp count_goldbach.cpp

int main()
{
	size_t MaxNumberOfThreads = 128;
	omp_set_num_threads(MaxNumberOfThreads);
	cout << "Number of Threads: " << MaxNumberOfThreads << endl;
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
    //double t1, tstart; //time measurement
    
	int N = (n-4)/2 + 1;        //Only consider even numbers from 4-n
    //vector<int> pair_cont(N,0); // all counts are initialized as zero
    vector<int> pair_contTwo(N,0);
    auto primes = get_primes(n); // get all the primes until n
       
    //tstart  = omp_get_wtime();
    //count_goldbach(n, pair_cont, primes); 
    //t1 = omp_get_wtime() - tstart;
       

//#################################################################### 
//Check the number with maximum pairs
    //cout << "The even number has most prime pairs is" <<" "<< 2*(2 + distance(pair_cont.begin(), max_element(pair_cont.begin(),pair_cont.end())))<< endl;
	//cout << "It has " << *max_element(pair_cont.begin(),pair_cont.end()) << " prime pairs" << endl;

//#################################################################### 
// Golbach Second Version
    double t2, tstart2; //time measurement
    tstart2  = omp_get_wtime();
    count_goldbachTwo(n, pair_contTwo, primes);   
    t2 = omp_get_wtime() - tstart2; 
//#################################################################### 
//Check the number with maximum pairs
    //cout << "The even number has most prime pairs is" <<" "<< 2*(2 + distance(pair_cont.begin(), max_element(pair_cont.begin(),pair_cont.end())))<< endl;
	//cout << "It has " << *max_element(pair_cont.begin(),pair_cont.end()) << " prime pairs" << endl;
//####################################################################    
// Performance evaluation
    cout << endl;
    cout.precision(4);
    //cout << "Total time	First Version: " << t1 << endl;
    cout << "Total time	Second Version: " << t2 << endl;
// Calculations CheckUp
	if(n<= 100000)
	{
		vector<int> contCheck(N,0);
		count_goldbachCheck(n,contCheck,primes);
		bool result = true;
		// First Version CheckUp
		//for(int k=0;k<N;++k){ if(pair_cont[k] != contCheck[k]){result = false; cout<< "Result First Version: Incorrect!!" << endl; break;} }
		//if(result == true){ cout<< "Result Fist Version: Correct!" << endl; }	
		// Second version CheckUp
	    result = true;
	    for(int k=0;k<N;++k){ if(pair_contTwo[k] != contCheck[k]){result = false; cout<< k <<  "Result Second Version: Incorrect!!" << endl;break;} }
		if(result == true){ cout<< "Result Second Version: Correct!" << endl; }
	}
    return 0;
}
