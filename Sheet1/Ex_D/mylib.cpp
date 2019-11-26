#include "mylib.h"
#include <cmath>
#include <cassert>       // assert()
#include <vector>
#include <iostream>
using namespace std;

double scalar(long long int n)
{
    double sum = 0.0;
    double eachsum = 0.0;
    for (int i = 1; i <= n; ++i)
    {
        eachsum = (1.0/i)*(1.0/i);
        sum += eachsum;

    }
    return sum;
}



