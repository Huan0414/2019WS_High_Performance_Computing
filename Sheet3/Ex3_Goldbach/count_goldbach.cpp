#include "mylib.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>

using namespace std;

void count_goldbach(const int n, vector<int> &pair_cont, vector<int> const &primes)
{
    
    int primeN = primes.size();
    //cout << "primeN: " << primeN << endl;
	pair_cont[0] = 1; // 4 has an pair (2+2)
	int sum = 0;
	
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


}
