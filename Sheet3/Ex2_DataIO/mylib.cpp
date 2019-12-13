#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

void cal_min_max_mean(vector<double> const &a, vector<double> &b)
{
 
 //Initialization
    size_t const N = a.size();
    double dAritSum = 0.0;
    double dGeoSum = 1.0;
    double dHarSum = 0.0;

    double dAritMean;
    double dGeoMean;
    double dHarMean;
      
// determine minimum and maximum
    b.at(0) = *min_element(a.begin(),a.end()); //minimum
    b.at(1) = *max_element(a.begin(),a.end()); //maximum 

// Calculate 'sum' value
    for (size_t i = 0; i < N; ++i)
    {
        dAritSum += a[i];
        dGeoSum *= pow(a[i], 1.0/N); 
        dHarSum += 1/a[i];
    }

// Calculate Mean Values
    b.at(2) = dAritSum/N; 			 //arithmetic mean
    b.at(3)  = dGeoSum;				 //geometric mean
    b.at(4)  = N/dHarSum;			 //harmonic mean
 
// Calculate standard deviation
    double dVarSum = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
    dVarSum += (a[i] - b.at(2))*(a[i] - b.at(2));
    }
    b.at(5) = sqrt(dVarSum/N); //standard deviation

}





