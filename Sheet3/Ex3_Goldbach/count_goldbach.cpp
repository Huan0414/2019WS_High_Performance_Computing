#include "mylib.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <math.h> 

using namespace std;

void count_goldbachCheck(const int n, vector<int> &pair_cont, vector<int> const &primes)
{
    int primeN = primes.size();
	pair_cont[0] = 1; // 4 has a pair (2+2)
	int sum = 0;
    
    for (int i=1; i< primeN; ++i)
	{		
		for (int j=i; j < primeN; ++j)
		{
			sum = primes[i] + primes[j];
			if(sum > n){break;}
			else {pair_cont[sum/2-2] += 1;}
		}
    }
    
    //cout << "checkUp version result: " ;
    //for(int k=0; k<(n-4)/2 + 1; ++k){cout<< pair_cont[k]<< ", ";}
    //cout << endl;
}

void count_goldbachTwo(const int n, vector<int> &pair_cont, vector<int> const &primes)
{
	int N = (n-4)/2 + 1;        //Only consider even numbers from 4-n
    int primeN = primes.size();
	int sum = 0;
	float const PercentageInnerForLoop = 0.8;

	double upperBound = primeN-1;
	vector<int> upperBounds(primeN,primeN-1);
	
	// Calculation of bounds
	
	int indexOfPrimeLowerEqualThanEnHalf;
		
	//double t1, tstart; //time measurement	for boundaries calculation
	//tstart  = omp_get_wtime();
	// Calculation of Boundaries
	for (int i=1; i< primeN; ++i) 
	{	
		if(primes[i] <= n/2)
		{
			for (int j=upperBounds[i-1]; j>=i; --j)
			{
				//cout << "i, j: " << i << " " << j << "primes[i], primes[j]: " << primes[i] << " " << primes[j] << endl;
				sum = primes[i] + primes[j];
				//cout << "sum: " << sum << endl;
				if(sum <= n){upperBounds[i+1] = upperBounds[i]; break;}
				else{upperBounds[i] -= 1;}
			}
		}else{ indexOfPrimeLowerEqualThanEnHalf = i-1; break; }
    }
    
    //for (int i=1; i< primeN; ++i){cout<< upperBoundsTest[i] << ",";}
    //cout << endl;
    //t1 = omp_get_wtime() - tstart;
    //cout << "Total time	Boundaries: " << t1 << endl;
    
    // Processing of 100% - PercentageInnerForLoop of the data with parallelization of the outer loop
    #pragma omp parallel for default(none) shared(cout,primes,primeN,n,upperBounds,indexOfPrimeLowerEqualThanEnHalf,PercentageInnerForLoop) private(sum) reduction(VecAdd:pair_cont) schedule(dynamic,10) 
	for (int i= (int)indexOfPrimeLowerEqualThanEnHalf*PercentageInnerForLoop + 1; i <= indexOfPrimeLowerEqualThanEnHalf; ++i)
	{	
		for (int j=i; j <= upperBounds[i]; ++j)
		{
			sum = primes[i] + primes[j];
			//assert(sum/2-2 < N);
			pair_cont[sum/2-2] += 1;
		}
    }
    
    // Processing PercentageInnerForLoop of the data with paralelization of the inner for loop  
	for (int i=1; i <= indexOfPrimeLowerEqualThanEnHalf*PercentageInnerForLoop; ++i)
	{	
		#pragma omp parallel for default(none) shared(primes,n,pair_cont,i,upperBounds) private(sum)
		for (int j=i; j <= upperBounds[i]; ++j)
		{
			sum = primes[i] + primes[j];
			//if(sum > n){cout<< "!!! sum > n ----- i,j: " << i << "," << j << " -- primes[i], primes[j]: " << primes[i] << "," << primes[j] << endl;}
			//assert(sum/2-2 < N);
			pair_cont[sum/2-2] += 1;
		}
    }
    
    pair_cont[0] = 1; // 4 has an pair (2+2)
    
}

void count_goldbach(const int n, vector<int> &pair_cont, vector<int> const &primes)
{  
    int N = (n-4)/2 + 1;        //Only consider even numbers from 4-n
    int primeN = primes.size();
    //cout << "primeN: " << primeN << endl;
	vector<int> contCheck(N,0);
	//contCheck[0] = 1; // 4 has a pair (2+2)
	int sum = 0;
	double partition = 0;
    
    #pragma omp parallel for default(none) shared(primes,primeN,n) private(sum) reduction(VecAdd:pair_cont) schedule(dynamic,10)
	for (int i=1; i< primeN; ++i)
	{		
		for (int j=i; j < primeN; ++j)
		{
			sum = primes[i] + primes[j];
			if(sum > n){break;}
			else {pair_cont[sum/2-2] += 1;}
		}
    }
    
    pair_cont[0] = 1; // 4 has an pair (2+2), dont move to the beginning! if the code is parallelized, each thread will add one to the pairs of 4
    
    //cout << "first version result: " ;
    //for(int k=0; k<N; ++k){cout<< pair_cont[k]<< ", ";}
    //cout << endl;
    
	//// Partition of sections into equal-computation segments (no real difference of performance with dynamic 10 old code)
	//cout << "New code" << endl;
	//int PartitionSize = (int) n/1000.0;
	////assert(PartitionSize>9);
	//if( PartitionSize < 10 ){PartitionSize = 10;}
	//int NumberOfThreads = 64;
	//cout << "PartitionSize: " << PartitionSize << endl;
	//cout << "Number of Threads: " << NumberOfThreads << endl;
	//vector<int> partitionOfPrimesForThreads(PartitionSize+1);
	//partitionOfPrimesForThreads[0] = 1;
	//partitionOfPrimesForThreads[PartitionSize] = primeN + 1;
	
	//for(int k=1; k<PartitionSize; ++k)
	//{	partition = floor( primeN*(1.0-sqrt(( PartitionSize*1.0-k )/( PartitionSize*1.0 ))) );
		//partitionOfPrimesForThreads[k]= (int) partition;}
		////cout<<"print value of x_k: " << partitionOfPrimesForThreads[k] << endl;
	 
	//omp_set_num_threads(NumberOfThreads);
	//#pragma omp parallel for default(none) shared(primes,primeN,n,PartitionSize,partitionOfPrimesForThreads) private(sum) reduction(VecAdd:pair_cont) schedule(dynamic)
	//for(int thread = 0; thread<PartitionSize; ++thread)
	//{
		//for (int i=partitionOfPrimesForThreads[thread]; i< partitionOfPrimesForThreads[thread+1]; ++i)
		//{		
				//for (int j=i; j < primeN; ++j)
				//{
					//sum = primes[i] + primes[j];
					//if(sum > n){break;}
					//else {pair_cont[sum/2-2] += 1;}
				//}
		//}
    //}
    
}
