#include "mylib.h"
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

//int main(int argc, char** argv)
int main()
{
    int const NROW = 4, MCOL = 3;           // initialize constants

    vector<double> const u({1, 2, 3});      // initialize u

// ---------------------------------------------------------------------
// matrix as 2D vector
    vector<vector<double>> const M({ {4, -1, -0.5}, { -1, 4, -1}, { -0.5, -1, 4}, {3, 0, -1} });

    assert( NROW ==  static_cast<int>(M.size()) );        // check number of rows
    assert( MCOL ==  static_cast<int>(M.back().size()) ); // check number of columns in last row

    //auto f1 = MatVec(M,u);
    auto f1 = M*u;

    for (size_t k=0; k<f1.size(); ++k)
    {
        cout << f1[k] << "  ";
    }
    cout << endl;
// -----------------------------------------------------------------------------
// matrix as 1D vector
    vector<double> const S({4, -1, -0.5, -1, 4, -1, -0.5, -1, 4, 3, 0, -1 });
    assert( NROW * MCOL == static_cast<int>(S.size()) );

    //auto f2 = MatVec(S,u);
    auto f2 = S*u;

    for (size_t k=0; k<f2.size(); ++k)
    {
        cout << f2[k] << "  ";
    }
    cout << endl;

  return 0;
}
