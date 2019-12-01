#include <vector>
#include <iostream>
using namespace std;

double scalar(vector<double> const &a, long long int const N)
{
    float sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
        sum += a[i];
    }
    return sum;
}

double Kahan_scalar(vector<double> const &a, long long int const N)
{
    // accumulator
    float sum = 0.0;
    // correction
    float c = 0.0;
    float y;
    float t;

    for (int i = 0; i < N; ++i)
    {
        y = a[i] - c;
        //
        t = sum + y;
        // new correction with lower part of y, will be added next loop
        c = (t - sum) - y;
        // non-accurate summation
        sum = t;
        
        //t = sum + a[i];
        //if (sum >= a[i]) {c += (sum- t) + a[i];}
        //else {c += (a[i] - t) + sum;}
        //sum = t; 
    }
    
    return sum+c;
}


