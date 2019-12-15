#include <cmath>
#include <vector>
#include <algorithm>
#include "mylib.h"
#include <omp.h>
using namespace std;

//void cal_min_max_mean(vector<double> const &a), vector<double> &b)
vector<double> cal_min_max_mean(vector<double> const &a)
{
 //Initialization
    size_t const N = a.size();
	vector<double> b(6);
    b[0] = a[0];
    b[1] = a[0];
    b[2] = 0.0;
    b[3] = 1.0;
    b[4] = 0.0;
    b[5] = 0.0;
    //double dAritSum = 0.0;
    //double dGeoSum = 1.0;
    //double dHarSum = 0.0;
    //double dVarCal  = 0.0;
    //double dAritMean,dGeoMean,dHarMean,dVarSum;    
      
// determine minimum and maximum
    //b.at(0) = *min_element(a.begin(),a.end()); //minimum
    //b.at(1) = *max_element(a.begin(),a.end()); //maximum

// Calculate 'sum' value
	#pragma omp parallel for default(none) shared(a,N) reduction(ReductionForMeansCalculation:b)
    
    for (size_t i = 0; i < N; ++i)
    {
        if(b[0] > a[i]){b[0] = a[i];}
        if(b[1] <a[i]){b[1] = a[i];}
        b[2] += a[i];
        b[3] *= pow(a[i], 1.0/N); 
        b[4] += 1.0/a[i];
        b[5] += a[i]*a[i];
        //dAritSum += a[i];
        //dVarCal +=  a[i]*a[i];
        //dGeoSum *= pow(a[i], 1.0/N); 
        //dHarSum += 1.0/a[i];      
    }
	
	// Calculate Mean Values
    b[2] /= N; 			 //arithmetic mean
    //b[3] = dGeoSum;				 //geometric mean
    b[4]  = N/b[4];			 //harmonic mean
// Calculate standard deviation
	b[5] = b[5] - N*b[2]*b[2];
	b[5] = sqrt(b[5]/(N-1)); //standard deviation 
	return b;
// Calculate Mean Values
    //b.at(2) = dAritSum/N; 			 //arithmetic mean
    //b.at(3)  = dGeoSum;				 //geometric mean
    //b.at(4)  = N/dHarSum;			 //harmonic mean
// Calculate standard deviation
	//dVarSum = dVarCal - 2.0*dAritSum*b.at(2) + N*b.at(2)*b.at(2);
	//b.at(5) = sqrt(dVarSum/(N-1)); //standard deviation 
//// Calculate standard deviation
    //double dVarSum = 0.0;
    //for (size_t i = 0; i < N; ++i)
    //{
    //dVarSum += (a[i] - b.at(2))*(a[i] - b.at(2));
    //}
}





