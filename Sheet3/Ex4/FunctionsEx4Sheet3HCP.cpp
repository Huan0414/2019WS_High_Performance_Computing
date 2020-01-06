#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include "mylib.h"
#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <sstream>
#include <string>
#include <assert.h>
#include <cmath>
#include <vector>

using namespace std;

inline
double InnerProduct(vector<double> const &x, vector<double> const &y, double &sum)
{
    assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    //double 
    sum = 0.0;
    //#pragma omp parallel for default(none) shared(x,y,N) reduction(+:sum)
    for (size_t i = 0; i < N; ++i)
    {
        //int th_id  = omp_get_thread_num();
        sum += x[i] * y[i];
        //sum += exp(x[i])*log(y[i]);
    }
    return sum;
}

void MatrixVectorMultiplication(vector<vector<double>> const &A, vector<double> const &x, vector<double> &product)
{
    int const mrows = A.size();         // #matrix rows

    #pragma omp parallel for default(none) shared(A,x,product,mrows)
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
  //#pragma omp parallel for default(none) shared(A,Transpose,RowsA,ColumnsA) collapse(2)
  for(int j=0; j<ColumnsA; ++j){
    for(int i=0; i<RowsA; ++i){
      Transpose[j][i] = A[i][j];
    }
  }
}

void MatrixMultiplication( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   assert(columnA == rowB); // Checks if dimensions match
   std::vector<std::vector<double>> TransposeB(columnB,std::vector<double>(rowB));
   TransposeMatrix(b,TransposeB);
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   //#pragma omp parallel for default(none) shared(a,product,rowA,columnB,TransposeB)
   for(int i=0; i<rowA; ++i){
     for(int j=0; j<columnB; ++j){
       InnerProduct(a[i],TransposeB[j],product[i][j]);
       //product[i][j] = 0.0;
       //for(int k=0; k<columnA; ++k){
       //product[i][j]+=a[i][k]*b[k][j];}
   }}
}

void CheckEquality(std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> &B){
  int RowsA=A.size(), ColumnsA=A[0].size(), RowsB=B.size(), ColumnsB=B[0].size();
  assert(RowsA == ColumnsB);                  // Dimension Check
  assert(ColumnsA == RowsB);            // Dimension Check
  for(int i=0; i<RowsA; ++i){
    for(int j=0; j<ColumnsA; ++j){
      assert(fabs(B[j][i]- A[i][j]) < 1e-3 ); }}
}

// Assumes the inizialization of the product matrix as a zero matrix.
void MatrixMultiplicationParallel( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   double variable1, variable2;
   assert(columnA == rowB); // Checks if dimensions match
   //std::vector<std::vector<double>> TransposeB(columnB,std::vector<double>(rowB));
   //TransposeMatrix(b,TransposeB);
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   #pragma omp parallel for default(none) shared(a,b,product,rowA,columnB,columnA) private(variable1,variable2) //collapse(2)
   for(int i=0; i<rowA; ++i){ 
     for(int k=0; k<columnA; ++k){
		//variable1 = a[i][k];
		//#pragma omp parallel for default(none) shared(a,b,i,k,variable1,product,columnB) //collapse(2)
		for(int j=0; j<columnB; ++j){
			//#pragma omp atomic read
			//variable1 = a[i][k];
			//#pragma omp atomic read
			//variable2 = b[k][j];
			//#pragma omp atomic write
			//product[i][j]+= variable1*variable2;
			product[i][j]+= a[i][k]*b[k][j];
   }}}
   
}

void PolynomialEvaluation(std::vector<double> const &a, std::vector<double> const &x, std::vector<double> &result){
    int Xsize = x.size();
    int Asize = a.size();
    int Rsize = result.size();
    assert(Rsize == Xsize);

    #pragma omp parallel for default(none) shared(a,x,Xsize,Asize,result)
    for (int j = 0; j<Xsize;j++)
    {
        result[j]=0.0;
        double x_iterative = 1.0;
        result[j] += a[0]*x_iterative;  // i = 0
        for (int i = 1; i<Asize;i++)
        {
            x_iterative *= x[j]; result[j] += a[i]*x_iterative;
        }
    }

}
