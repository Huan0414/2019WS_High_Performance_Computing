#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <ctime>
#include "mylib.h"
using namespace std;

int main(int argc, char **argv)
{
    int const NLOOPS = 50;        // chose a value such that the benchmark runs at least 10 sec.
    unsigned int N = 5000001;

//  Memory allocation
    cout << "Memory allocation\n";

    vector<double> x(N);

    cout.precision(2);
    cout << 2.0 * N *sizeof(x[0]) / 1024 / 1024 / 1024 << " GByte Memory allocated\n";
    cout.precision(6);
    cout << endl;

//##########################################################################
//  Input Parameters
//  case 1: input (1,4,16)
//  case 2: input (2,3,5)
//  case 3: input (1000,4000,16000)
//  case 4: With STL vector of arbitrary length containing the input data
    cout << "This exercise is to calculate arithmetic mean value, geometric "
         << "mean value and harmonic mean value -- By Huan CHEN\n";
    cout << endl;
    cout << "Case Explanation:\n";
    cout << "case 1: input (1,4,16).\n";
    cout << "case 2: input (2,3,5).\n";
    cout << "case 3: input (1000,4000,16000).\n";
    cout << "case 4: With STL vector of arbitrary length containing the input data.\n";
    cout << endl;

    int nCaseNr;

// Select case number
    cout << "Please Enter case number:";
    cin  >> nCaseNr;
    cout << "\nStart Benchmarking\n"<<endl;
    double tstart, t1;  // timer
// Do calculation
    tstart = clock();       // start timer

// Define input parameter array

    switch (nCaseNr)
    {
    case 1: x = {1,4,16};
            break;
    case 2: x = {2,3,5};
            break;
    case 3: x = {1000,4000,16000};
            break;
    case 4: for(int nNrValues = 0;nNrValues < 50000;nNrValues++)
            {
            double xValue;
            cout << "Please enter an arbitrary negative value when vector values are all entered"<<endl;
            cout << "Enter next number:"<<endl;
            cin >> xValue;
            if (xValue>0)
            {
                x[nNrValues] = xValue;
            }
            break;
    default: cout << "Please check the case number again!!"<<endl;
            break;
    }

// Print the result
    vector<double> dMeanArray(3);
    dMeanArray = calMean(x);
    cout << "The arithmetic mean value is:" << dMeanArray[0] <<endl;
    cout << "The geometric mean value is:" << dMeanArray[1] <<endl;
    cout << "The harmonic mean value is:" << dMeanArray[2] <<endl;

    t1 = clock() - tstart;  // stop timer
    t1 /= CLOCKS_PER_SEC;   // now, t1 in seconds
    t1 /= NLOOPS;           // divide by number of function calls

    return 0;
}
