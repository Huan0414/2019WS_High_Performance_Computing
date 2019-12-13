#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
vector<long> get_primes(long max);

int single_goldback(long n)
{

    //initialize counting number
      int deNr = 0;
    //get all the primes
    vector<long> primes = get_primes(n);
    //Define first prime
    vector<long>::iterator it;
    for(it = primes.begin(); it <= primes.end(); it++)
    {
        // Define the second prime number according to the first one
        long secPrime = n - *it;
        // Check if the second prime is within the primes vector
        if (count(primes.begin(),primes.end(),secPrime))
        {
            deNr++;
            // for condition like 10 = 5 + 5.
            if (secPrime == *it && count(primes.begin(),primes.end(),secPrime))
            {
                deNr ++;
            }
        }
    }
    deNr = deNr/2;
    return deNr;
}
