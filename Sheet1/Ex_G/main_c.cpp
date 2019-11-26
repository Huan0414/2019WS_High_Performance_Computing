#include <cassert>
#include <iostream>
#include <vector>
#include <ctime>
#include "DenseMatrix.h"
using namespace std;

int main()
{

    DenseMatrix const syMatrix(1000,1000);
    vector<double> const u{{1,2,3}};

    int const NLOOPS=100; // the overall code should run approx. 10 sec.

    //Mult Function
    vector<double> f1;
    double t1 = clock(); // start timer
    for (int k=1; k<NLOOPS; ++k) {f1 = syMatrix.Mult(u);}
    t1 = (clock()-t1)/CLOCKS_PER_SEC/NLOOPS;

    //MultT Function
    vector<double> f2;
    double t2 = clock(); // start timer
    for (int k=1; k<NLOOPS; ++k) {f2 = syMatrix.MultT(u);}
    t2 = (clock()-t2)/CLOCKS_PER_SEC/NLOOPS;

    cout<<"Time for one loop Mult function is:"<<t1 <<" "<<"seconds"<<endl;
    cout<<"Time for one loop MultT function is:"<<t2 <<" "<<"seconds"<<endl;

    return 0;
}
