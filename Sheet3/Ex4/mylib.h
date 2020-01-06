#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <cmath>
#include <vector>

using namespace std;

// g++ -Ofast -fopenmp Ex4Sheet3.cc && ./a.out

// Function InnerProduct
double InnerProduct(vector<double> const &x, vector<double> const &y, double &sum);

// Function Matrix vector multiplication
void MatrixVectorMultiplication(vector<vector<double>> const &A, vector<double> const &x, vector<double> &product);


// Function TransposeMatrix
void TransposeMatrix(std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> &Transpose);

// Function Matrix Multiplication, not parallelized. 
void MatrixMultiplication( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product);
                          
// Function CheckEquality: checks if two matrices are equal entry by entry up to a difference of 10 minus 3, 
// an assert jumps if the difference goes beyond this threshold
void CheckEquality(std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> &B);

// Parallelized version of Matrix Multiplication. Assumes the inizialization of the product matrix as a zero matrix. 
void MatrixMultiplicationParallel( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product); 
                          
// Function Polynomial evaluation with parallelization
void PolynomialEvaluation(std::vector<double> const &a, std::vector<double> const &x, std::vector<double> &result);
