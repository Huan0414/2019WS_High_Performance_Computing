#include "mylib.h"
#include <cassert>
#include <vector>
using namespace std;


vector<double> MatVec(vector<vector<double>> const &A, vector<double> const &u)
{
    int const nrows = static_cast<int>(A.size());         // #matrix rows
    int const mcols = static_cast<int>(A[0].size());      // #matrix columns
    assert( mcols ==  static_cast<int>(u.size()) );       // check compatibility inner dimensions

    vector<double> f(nrows);                 // allocate resulting vector

                                             // matrix times vector: f := A * u
    for (int i = 0; i < nrows; ++i) {
        double tmp = 0.0;			         // initialize f[i]
        for (int j = 0; j < mcols; ++j) {
            tmp = tmp + A[i][j] * u[j];
        }
        f[i] = tmp;
        //cout << A[i].data() << endl;           // Address of a[i][0]
    }

    return f;
}
// ---------------------------------------------------------------------

vector<double> MatVec(vector<double> const &A, vector<double> const &u)
{
    int const nelem = static_cast<int>(A.size());      // #elements in matrix
    int const mcols = static_cast<int>(u.size());      // #elements in vector <==> #columns in matrix

    assert(nelem % mcols == 0);                        // nelem has to be a multiple of mcols (==> #rows)
    int const nrows = nelem/mcols;                     // integer division!

    vector<double> f(nrows);                 // allocate resulting vector
                                             // matrix times vector: f := A * u
    for (int i = 0; i < nrows; ++i) {
        double tmp = 0.0;			         // initialize f[i]
        for (int j = 0; j < mcols; ++j) {
            tmp = tmp + A[i*mcols+j] * u[j];
        }
        f[i] = tmp;
        //cout << A[i*mcols].data() << endl;           // Address of a[i][0]
    }

    return f;
}
