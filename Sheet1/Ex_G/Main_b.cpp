#include <iostream>
#include <cassert>
#include <vector>
#include <math.h>
#include "DenseMatrix.h"
using namespace std;

//***************************************************************
//main function
int main()
{
DenseMatrix const M(5,3); // Dense matrix, also initialized

vector<double> const u{{1,2,3}};
vector<double> f1 = M.Mult(u);

vector<double> const v{{-1,2,-3,4,-5}};
vector<double> f2 = M.MultT(v);

return 0;
}
