#ifndef FUNCTIONSEX3SHEET2_H_INCLUDED
#define FUNCTIONSEX3SHEET2_H_INCLUDED

#include<vector>
#include <limits>
#include <cmath>
#define M_PI    3.14159265358979323846 /*pi*/

// Function Inner Product
void InnerProduct(std::vector<double> const &vecOne, std::vector<double> const &vecTwo, double &product);

// L2 Norm of vector
void LTwoNorm(std::vector<double> const &vec, double &product);

// Function Matrix-Vector Multiplication
void MatrixVectorMultiplication(std::vector<std::vector<double>> const &A, std::vector<double> const &x, std::vector<double> &product);

// Function TransposeMatrix
void TransposeMatrix(std::vector<std::vector<double>> const &A, std::vector<std::vector<double>> &Transpose);

// Function MatrixMultiplication
void MatrixMultiplication( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b, std::vector<std::vector<double>> &product);

// Function MatrixMultilplicationTranspose , optimization of matrix Multiplication by transposing
void MatrixMultiplicationTranspose( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b, std::vector<std::vector<double>> &product);

// Function MatrixMultiplication Column-wise access
void MatrixMultiplicationColumnwiseAccess( std::vector<std::vector<double>> const &a, std::vector<std::vector<double>> const &b, std::vector<std::vector<double>> &product);

// Function Polynomial Evaluation
void PolynomialEvaluation(std::vector<double> const &a, std::vector<double> const &x, std::vector<double> &result);

double scalar(std::vector<double> const &a, long long int const N);

double Kahan_scalar(std::vector<double> const &a, long long int const N);

#endif // FUNCTIONSEX3SHEET2_H_INCLUDED
