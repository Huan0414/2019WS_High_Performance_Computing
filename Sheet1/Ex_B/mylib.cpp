#include <iostream>
#include "mylib.h"
#include <cmath>
#include <cassert>       // assert()
#include <vector>
using namespace std;

vector<double> calMean(vector<double> const &x)
{
    size_t const N = x.size();

    //Initialization
    double dAritSum = 0.0;
    double dGeoSum = 1.0;
    double dHarSum = 0.0;

    double dAritMean;
    double dGeoMean;
    double dHarMean;

// Calculate 'sum' value
    for (size_t i = 0; i < N; ++i)
    {
        dAritSum += x[i];
        dGeoSum *= x[i];
        dHarSum += 1/x[i];
    }

// Calculate Mean Values
    dAritMean = dAritSum/N;
    dGeoMean  = pow(dGeoSum,1.0/N);
    dHarMean  = N/dHarSum;

// Calculate standard deviation
    double dVarSum = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
    dVarSum += (x[i]-dAritMean)*(x[i]-dAritMean);
    }
    double dSigama = sqrt(dVarSum/N);

// Store mean values into array and return to main function
    vector<double> dMeanArray = {dAritMean,dGeoMean,dHarMean,dSigama};
    return  dMeanArray;
}


