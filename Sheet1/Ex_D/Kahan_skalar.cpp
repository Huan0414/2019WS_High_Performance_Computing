#include "mylib.h"
#include <cmath>
#include <cassert>       // assert()
#include <vector>
using namespace std;

double Kahan_scalar(long long int n)
{
    // accumulator
    double sum = 0.0;
    // correction
    double c = 0.0;
    double y;
    double t;

    for (int i = 1; i <= n; ++i)
    {
        //
        y = (1.0/i)*(1.0/i) - c;
        //
        t = sum + y;
        // new correction with lower part of y, will be added next loop
        c = (t - sum) - y;
        // non-accurate summation
        sum = t;
    }

    return sum;
}
