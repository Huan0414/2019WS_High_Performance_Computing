#include "functionsEx3Sheet2.h"

#include<assert.h>
#include<cmath>
#include<vector>

// Function Inner Product
void InnerProduct(std::vector<double> const &vecOne, std::vector<double> const &vecTwo, double &product)
{
	assert(vecOne.size() == vecTwo.size());
    product = 0.0;
    for (int i = 0; (unsigned)i <= vecOne.size()-1; i++){
      product = product + (vecOne[i])*(vecTwo[i]);
    }
    //return product;
}

// L2 Norm of vector
void LTwoNorm(std::vector<double> const &vec, double &product)
{
    InnerProduct(vec,vec,product);
    // product = 0.0;
    //for (int i = 0; (unsigned)i <= vec.size()-1; i++){
    //  product = product + (vec[i])*(vec[i]);
    //}
    product = sqrt(product);
    //return product;
}

// Function Matrix-Vector Multiplication
void MatrixVectorMultiplication(std::vector<std::vector<double>> const &A, std::vector<double> const &x, std::vector<double> &product)
{
    int const mrows = A.size();         // #matrix rows

    for (int i = 0; i < (unsigned)mrows; ++i) {
      InnerProduct(A[i],x,product[i]);                    // inner product of row i with vector x
      //cout<< product[i] << " ";
    }
}

// Function TransposeMatrix
void TransposeMatrix(std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> &Transpose){
  int RowsA=A.size(), ColumnsA=A[0].size(), RowsT=Transpose.size(), ColumnsT=Transpose[0].size();
  assert(RowsA == ColumnsT);                  // Dimension Check
  assert(ColumnsA == RowsT);            // Dimension Check
  for(int i=0; i<RowsA; ++i){
    for(int j=0; j<ColumnsA; ++j){
      Transpose[j][i] = A[i][j];
    }
  }
}

// Function MatrixMultiplication
void MatrixMultiplication( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   assert(columnA == rowB); // Checks if dimensions match
   //vector<vector<double>> TransposeB(columnB,vector<double>(rowB));
   //TransposeMatrix(b,TransposeB);
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   for(int i=0; i<rowA; ++i){
     for(int j=0; j<columnB; ++j){
       //InnerProduct(a[i],TransposeB[j],product[i][j]);
       product[i][j] = 0.0;
       for(int k=0; k<columnA; ++k){
       product[i][j]+=a[i][k]*b[k][j];}
   }}
}

void MatrixMultiplicationTranspose( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   assert(columnA == rowB); // Checks if dimensions match
   std::vector<std::vector<double>> TransposeB(columnB,std::vector<double>(rowB));
   TransposeMatrix(b,TransposeB);
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   for(int i=0; i<rowA; ++i){
     for(int j=0; j<columnB; ++j){
       InnerProduct(a[i],TransposeB[j],product[i][j]);
       //product[i][j] = 0.0;
       //for(int k=0; k<columnA; ++k){
       //product[i][j]+=a[i][k]*b[k][j];}
   }}
}

void PolynomialEvaluation(std::vector<double> const &a, std::vector<double> const &x, std::vector<double> &result){
    int Xsize = x.size();
    int Asize = a.size();
    int Rsize = result.size();
    assert(Rsize == Xsize);

    for (int j = 0; j<Xsize;j++)
    {
        result[j]=0.0;
        double x_iterative = 1.0;
        for (int i = 0; i<Asize;i++)
        {
            if (i==0) {result[j] += a[i]*x_iterative;}
            else {x_iterative *= x[j]; result[j] += a[i]*x_iterative;}
        }
    }

}

double scalar(std::vector<double> const &a, long long int const N)
{
    float sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
        sum += a[i];
    }
    return sum;
}

double Kahan_scalar(std::vector<double> const &a, long long int const N)
{
    // accumulator
    float sum = 0.0;
    // correction
    float c = 0.0;
    float y;
    float t;

    for (int i = 0; i < N; ++i)
    {
        y = a[i] - c;
        //
        t = sum + y;
        // new correction with lower part of y, will be added next loop
        c = (t - sum) - y;
        // non-accurate summation
        sum = t;

        //t = sum + a[i];
        //if (sum >= a[i]) {c += (sum- t) + a[i];}
        //else {c += (a[i] - t) + sum;}
        //sum = t;
    }

    return sum+c;
}
