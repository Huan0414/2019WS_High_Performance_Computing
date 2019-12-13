#include <vector>
using namespace std;

#include <mayer_primes.h>

void count_goldback(const int n, vector<int> &pair_cont)
{
    auto primes = get_primes(n); // get all the primes until n
    int primeN = primes.size();

	int k; 
	pair_cont[0] = 1; // 4 has an pair (2+2)
    for (int i=1; i< primeN ; ++i)
    {		
			for (int j=i; j < primeN; j++)
			{
				int sum = primes[i] + primes[j];
				if (sum>=4 && sum<=n) {pair_cont[sum/2-2] += 1;}
			}
    }

}
