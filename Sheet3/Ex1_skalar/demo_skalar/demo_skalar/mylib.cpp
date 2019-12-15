#include "mylib.h"
#include <cassert>       // assert()
#include <iostream>
#include <numeric>       // iota()
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
using namespace std;

double scalar(vector<double> const &x, vector<double> const &y)
{
    assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    double sum = 0.0;
	#pragma omp parallel for default(none) shared(x,y,N) reduction(+:sum)
    for (size_t i = 0; i < N; ++i)
    {
        sum += x[i] * y[i];
        //sum += exp(x[i])*log(y[i]);
    }
    return sum;
}


double scalarNofor(vector<double> const &x, vector<double> const &y)
{	
	
    assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    double sum = 0.0;
    
	#pragma omp parallel default(none) shared(x,y,N) reduction(+:sum)
        {	
			int Thread = omp_get_thread_num();
			int NumberOfThreads = omp_get_num_threads();
			
			int i = Thread*N/NumberOfThreads;
			int j=0;
			
			while(j < N/NumberOfThreads) // averagely split addition to all threads
			{
				sum += x[i+j] * y[i+j]; 
				j += 1;
			}				
			
			if(Thread == 1)		// dealing with leftover part of addition,e.g. N=500001, just 1 left
			{ 
				j = 0;
				while(N-N%NumberOfThreads+j<N)
				{
					sum += x[N-N%NumberOfThreads+j] * y[N-N%NumberOfThreads+j]; 
					j+=1;
				}
			}
		}  
			
	return sum;
}


vector<int> reduction_vec_append(int n)
{ 
    vector<int> vec;
    
	#pragma omp declare reduction (OurAppend : vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
	#pragma omp parallel default(none) shared(cout, n) reduction(OurAppend:vec)
    {
		vec.resize(n);
        int threadNumber = omp_get_thread_num();
        iota( vec.begin(),vec.end(), threadNumber);
        #pragma omp critical
		{
			cout << "print out vec of each thread: " << threadNumber << endl;
			cout << vec << endl;
		}
		
		//for (int t = 0; t < omp_get_num_threads(); t++) {
        //#pragma omp barrier
			//if (t == omp_get_thread_num()) {
            //vec.resize(n);
			//int threadNumber = omp_get_thread_num();
			//iota( vec.begin(),vec.end(), threadNumber);
			//#pragma omp critical
			//{cout << "print out vec of each thread: " << threadNumber << endl;
			//cout << vec << endl;}
			//}
		//}
        //#pragma omp critical
        //cout << omp_get_thread_num() << " : " << vec.size() << endl;
        
        
	}
    return vec;
}

vector<int> reduction_vec_append_manually(int n)
{ 
    //vector<int> vec;
    vector<int> global_vec;
	#pragma omp parallel default(none) shared(cout, n,global_vec) 
    {
        //vec.resize(n);
        vector<int> vec(n);
        int threadNumber = omp_get_thread_num();
        iota(vec.begin(),vec.end(), threadNumber);
        #pragma omp critical
		{
			cout << "print out vec of each thread: " << threadNumber << endl;
			cout << vec << endl;
		}
		
		// do appending manully
		for (int t = 0; t < omp_get_num_threads(); t++) {
            #pragma omp barrier
			if (t == omp_get_thread_num()) {
            global_vec.insert(global_vec.end(), vec.begin(), vec.end());
			}
		}
	}
    return global_vec;
}

vector<int> reduction_vec(int n)
{ 
    vector<int> vec(n);
#pragma omp parallel default(none) shared(cout) reduction(VecAdd:vec)
    {
        //#pragma omp critical
        //cout << omp_get_thread_num() << " : " << vec.size() << endl;
        iota( vec.begin(),vec.end(), omp_get_thread_num() );
    }
    return vec;
}
