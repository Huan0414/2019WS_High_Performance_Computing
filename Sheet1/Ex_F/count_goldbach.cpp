#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
vector<long> get_primes(long max);

vector<int> count_goldback( long n)
{
    vector<long> primes = get_primes(n); // get all the primes until n
    int primeN = primes.size();

    vector<int> cont(n-3,0); // all counts are initialized as zero
    for (int i=0; i< primeN ; ++i)
    {
        for (int j=i; j < primeN; ++j)
        {
          int  sum = primes[i] + primes[j];
            if(sum>=4&&sum<=n) {cont[sum-3] += 1;}
        }
    }
    return cont;
}
