#include "mylib.h"
#include <cassert>       // assert()
#include <iostream>
#include <numeric>       // iota()
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
using namespace std;

double scalar(vector<double> const &x, vector<double> const &y)
{
    assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    double sum = 0.0;
#pragma omp parallel for default(none) shared(x,y,N) reduction(+:sum)
    for (size_t i = 0; i < N; ++i)
    {
        sum += x[i] * y[i];
        //sum += exp(x[i])*log(y[i]);
    }
    return sum;
}




vector<int> reduction_vec(int n)
{ 
    vector<int> vec(n);
#pragma omp parallel default(none) shared(cout) reduction(VecAdd:vec)
    {
        //#pragma omp critical
        //cout << omp_get_thread_num() << " : " << vec.size() << endl;
        iota( vec.begin(),vec.end(), omp_get_thread_num() );
    }
    return vec;
}
