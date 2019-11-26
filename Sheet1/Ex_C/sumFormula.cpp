#include <iostream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

long long int sumFormula(int n)
{
    long long int sumFormula = 0;
     // Take floor of n/3, n/5, n/15
    long long int floThree = n/3;
    long long int floFive = n/5;
    long long int floFivete = n/15;

    // Calculate the sum
    sumFormula  = (3*floThree*(floThree + 1) + 5*floFive*(floFive + 1) - 15*floFivete*(floFivete + 1))/2;

    return sumFormula;

}


