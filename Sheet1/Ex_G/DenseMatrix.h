#ifndef DENSEMATRIX_H_INCLUDED
#define DENSEMATRIX_H_INCLUDED
#include <iostream>
#include <cassert>
#include <vector>
#include <math.h>
using namespace std;
//*********************************************************
//Define a class
class DenseMatrix
{
    int nrow, mcolum;
    vector<vector<double>> A;

    public:
    //build the constructor matrix
    DenseMatrix (int n,int m);
    // Mult A*u
    vector<double> Mult( vector<double> const &u)const
    {
        vector<double> f(nrow);                 // allocate resulting vector
        for (int i = 0; i < nrow; ++i) {
            double tmp = 0.0;			         // initialize f[i]
            for (int j = 0; j < mcolum; ++j) {tmp = tmp + A[i][j] * u[j];}
            f[i] = tmp;
        }
        return f;
    }

    //Mutt AT*v
    vector<double> MultT( vector<double> const &v)  const{
        vector<double> f(mcolum);                 // allocate resulting vector
        for (int i = 0; i < mcolum; ++i) {
            double tmp = 0.0;			         // initialize f[i]
            for (int j = 0; j < nrow; ++j) {tmp = tmp + A[j][i] * v[j];}
            f[i] = tmp;
        }
        return f;
    }
};

DenseMatrix::DenseMatrix(int n, int m)
: nrow(n),mcolum(m), A(nrow, vector<double>(mcolum))
{
    int nm = max(nrow,mcolum);
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<mcolum;j++)
        {
            double xi = (10.0*i)/(nm-1.0) - 5.0;
            double xj = (10.0*j)/(nm-1.0) - 5.0;
            double fxi = 1.0/(1.0 + exp(-xi));
            double fxj = 1.0/(1.0 + exp(-xj));
            A.at(i).at(j) = fxi * fxj;
        }
    }
}
#endif // DENSEMATRIX_H_INCLUDED
