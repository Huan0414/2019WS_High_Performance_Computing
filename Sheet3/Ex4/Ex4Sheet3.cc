#include<assert.h>
#include<cmath>
#include<vector>

void MatrixVectorMultiplication(std::vector<std::vector<double>> const &A, std::vector<double> const &x, std::vector<double> &product)
{
    int const mrows = A.size();         // #matrix rows

    for (int i = 0; i < (unsigned)mrows; ++i) {
      InnerProduct(A[i],x,product[i]);                    // inner product of row i with vector x
      //cout<< product[i] << " ";
    }
}
