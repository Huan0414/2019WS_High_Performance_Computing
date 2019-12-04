#include<iostream>
#include<vector>
#include<assert.h>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include "functionsEx3Sheet2.h"

using namespace std;

// Main
int main()
{
  double const DoubleSize = sizeof(double); // Size of double
  double tstart, t1;  // timer
  // ###########################################################################################
  { /* */// Inner Product
  cout << "------" << endl;
  double const NLOOPSInnerProd = 200.0; // InnerProduct

  // Test data for InnerProduct
  int const VectorSize = 40000000;
  vector<double> a(VectorSize);
  vector<double> b(VectorSize);

  //Initialization of vectors "a" and "b" for InnerProduct
  for (int i=0;i<VectorSize;i++){
    a[i] = (i)%219 + 1;
    b[i] = 1.0/a[i];
  }

  double productInner;

  tstart = clock(); // start timer

  // Do Calculation
  for (int i = 0; i < NLOOPSInnerProd; ++i){
    InnerProduct(a,b,productInner);
  }
  // End Calculation

  //############################################# Print Timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;   // now, t1 in seconds
  cout << "InnerProduct "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSInnerProd;           // divide by number of function calls
  cout.precision(6);
  cout << "N (size of vector): " << VectorSize << endl;
  cout << "Number of Loops: " << NLOOPSInnerProd << endl;
  cout << "Timing in sec. per loop: " << t1 << endl;
  cout << "Memory allocated: " << 2.0*VectorSize*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*VectorSize/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 4.0*VectorSize*DoubleSize/t1/1024/1024/1024 << endl;
  //############################################# -------------
  }
  // ###########################################################################################
  { /* */// Matrix-Vector Multiplication
  double const NLOOPSMatrixVector = 150.0; // Matrix-Vector

  cout << "------" << endl;
  // Test data for Matrix-Vector multiplication.
  int const MROW=40000, NCOL=2000;           // initialize constants

  vector<double> x(NCOL);      // initialize u
  vector<vector<double>> A(MROW, vector<double>(NCOL));

  for (int i=0;i<MROW;++i)
  {
    for (int j=0;j<NCOL;++j)
      {
        A[i][j] = (i+j)%219 + 1;
        if(i==17) {x[j] = 1.0/A[17][j];}
        //cout << A[i][j] << "  ";
      }
     //cout << endl;
  }

  vector<double> z(MROW);
  tstart = clock(); // start timer
  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixVector; ++i){
     MatrixVectorMultiplication(A,x,z);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Matrix-Vector Multiplication "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixVector; // divide by number of function calls
  cout.precision(2);
  cout << "M times N matrix: " << "M = " << MROW << ", N = " << NCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixVector << endl;
  cout << "Timing in sec. per loop: : " << t1 << endl;
  cout << "Memory allocated (in GB): " << (MROW*NCOL + MROW + NCOL)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*MROW*NCOL/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< (4.0*MROW*NCOL+MROW)*DoubleSize/t1/1024/1024/1024 << endl;
  // ############################################ ----------------------------
  }
  // ###########################################################################################
  { /* */ // Matrix-Matrix Multiplication
  cout << "------" << endl;
  double const NLOOPSMatrixMatrix = 5.0; // Matrix-Matrix

  int const QROW= 1000, BROW = 1200, BCOL= 1000;
  vector<vector<double>> Q(QROW, vector<double>(BROW));
  //vector<vector<double>> QTranspose(QROW, vector<double>(BROW));
  vector<vector<double>> B(BROW, vector<double>(BCOL));


  for (int i=0;i<QROW;i++) // initialization of Q
  {
    for (int j=0;j<BROW;j++)
      {
        Q[i][j] = (i+j)%219 + 1;
      }
  }
  for (int i=0;i<BROW;i++) // Initialization of B
  {
    for (int j=0;j<BCOL;j++)
      {
        B[i][j] = (i+j)%219 + 1;
      }
  }

  //TransposeMatrix(Q,QTranspose);

  // ############################################ Calculation
  vector<vector<double>> product(QROW, vector<double>(BCOL));
  tstart = clock(); // start timer

  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixMatrix; ++i){
    MatrixMultiplication(Q,B,product);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Matrix-Matrix Multiplication "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixMatrix; // divide by number of function calls
  cout.precision(2);
  cout << "Matrix Dimensions (M times L and L times N): M = " << QROW << ", L = " << BROW << ", N = " << BCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixMatrix << endl;
  cout << "Timing in sec. per loop: " << t1 << endl;
  cout << "Memory allocated (in GB): " << (QROW*BROW + BROW*BCOL+QROW*BCOL)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << (2.0*QROW*BROW*BCOL)/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< (4.0*QROW*BROW*BCOL+QROW*BCOL)*DoubleSize/t1/1024/1024/1024 << endl;
  // ############################################ ----------------------------
  }
  // ###########################################################################################
  { /* */// Polynomial Evaluation
  cout << "------" << endl;
  double const NLOOPSPolynomialEvaluation = 10.0;
  int polynomialSize = 200000;
  int evaluations = 5000;
  vector<double> coefficients(polynomialSize);
  vector<double> x(evaluations);
  // test data
  for(int i=0; i<evaluations; i++){
    x[i] = (i)%100 + 1;
  }

  for(int i=0; i<polynomialSize; i++){
    coefficients[i] = (i)%219 + 1;
  }
  // Do Calculation
  vector<double> result(evaluations);
  tstart = clock(); // start timer
  for (int i=0; i < NLOOPSPolynomialEvaluation; ++i){
    PolynomialEvaluation(coefficients,x,result);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Polynomial Evaluation "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSPolynomialEvaluation; // divide by number of function calls
  cout.precision(2);
  cout << "Size of vector x: " << evaluations << endl;
  cout << "Polynomial size (N) : " << polynomialSize << endl;
  cout << "Number of Loops: " << NLOOPSPolynomialEvaluation << endl;
  cout << "Timing in sec. per loop: " << t1 << endl;
  cout << "Memory allocated (in GB): " << (2.0*evaluations+polynomialSize)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << (3.0*polynomialSize*evaluations)/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 7.0*evaluations*polynomialSize*DoubleSize/t1/1024/1024/1024 << endl;
  // ############################################ ----------------------------
  }
  // ###########################################################################################
  { /* */// L2 Norm
  cout << "------" << endl;
  double const NLOOPSL2 = 500.0; // L2Norm

  // Test data for L2Norm
  int const VectorSize = 40000000;
  vector<double> a(VectorSize);

  double Norm;

  tstart = clock(); // start timer

  // Do Calculation
  for (int i = 0; i < NLOOPSL2; ++i){
    LTwoNorm(a,Norm);
  }
  // End Calculation

  //############################################# Print Timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;   // now, t1 in seconds
  cout << "L2 Norm "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSL2;           // divide by number of function calls
  cout.precision(2);
  cout << "Vector Size (N): " << VectorSize << endl;
  cout << "Number of Loops: " << NLOOPSL2 << endl;
  cout << "Timing in sec. per loop: " << t1 << endl;
  cout << "Memory allocated: " << 2*VectorSize*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*VectorSize/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 3.0*VectorSize*DoubleSize/t1/1024/1024/1024 << endl;
   //############################################# -------------
   }
   // ###########################################################################################
  {  /* */// Kahan Summation vs Normal Summation
  cout << "------" << endl;
  int const NLOOPS = 50;        // chose a value such that the benchmark runs at least 10 sec.
  long long int N = 40000000;
  vector<double> a(N);

  // vector initialiation
  for (int i=1;i<=N;i++)
  {
    a[i] = (1.0/i)*(1.0/i);
  }

  double tstart, t1;  // timer

  // Do calculation
  tstart = clock();       // start timer

  for (int i=0;i<NLOOPS;i++)
  {
	Kahan_scalar(a,N);
  }
  t1 = clock() - tstart;

  //##########################################################################
  // Timings  and Performance
  t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Kahan Summation" << endl;
  cout.precision(6);
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPS; // divide by number of function calls
  cout.precision(6);
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Sum size: " << N << endl;
  cout << "Loops: " << NLOOPS << endl;
  cout << "Memory allocated (in GB): " << 2.0*N*sizeof(double)/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 4.0*N/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 12.0*N*sizeof(double)/t1/1024/1024/1024 << endl;
  //----------------------------------------------------------------------------
  cout << "Normal summation result is:" << scalar(a,N)<< endl;
  cout << "Kahan summation result is:"  << Kahan_scalar(a,N)<< endl;
  cout << "The difference: KahanSum - NormalSum = " << Kahan_scalar(a,N) - scalar(a,N)<<endl;
  cout << "The difference: pi^2/6 - NormalSum = " << pow(M_PI,2.0)/6.0 - scalar(a,N)<< endl;
  cout << "The difference: pi^2/6 - KahanSum = " << pow(M_PI,2.0)/6.0 - Kahan_scalar(a,N)<<endl;}
  // ###########################################################################################
  { /* */// MatrixMultiplicationTranspose
  // ###########################################################################################
  cout << "------" << endl;
  double const NLOOPSMatrixMatrixTranspose = 5.0; // Matrix-Matrix

  int const QROW= 1000, BROW = 1200, BCOL= 1000;
  vector<vector<double>> Q(QROW, vector<double>(BROW));
  //vector<vector<double>> QTranspose(QROW, vector<double>(BROW));
  vector<vector<double>> B(BROW, vector<double>(BCOL));


  for (int i=0;i<QROW;i++) // initialization of Q
  {
    for (int j=0;j<BROW;j++)
      {
        Q[i][j] = (i+j)%219 + 1;
      }
  }
  for (int i=0;i<BROW;i++) // Initialization of B
  {
    for (int j=0;j<BCOL;j++)
      {
        B[i][j] = (i+j)%219 + 1;
      }
  }

  //TransposeMatrix(Q,QTranspose);

  // ############################################ Calculation
  vector<vector<double>> product(QROW, vector<double>(BCOL));
  tstart = clock(); // start timer

  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixMatrixTranspose; ++i){
    MatrixMultiplicationTranspose(Q,B,product);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = clock() - tstart;
  t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Matrix-Matrix Multiplication Optimized with transposition "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixMatrixTranspose; // divide by number of function calls
  cout.precision(2);
  cout << "Matrix Dimensions (M times L and L times N): M = " << QROW << ", L = " << BROW << ", N = " << BCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixMatrixTranspose << endl;
  cout << "Timing in sec. per loop: " << t1 << endl;
  cout << "Memory allocated (in GB): " << (QROW*BROW + 2.0*BROW*BCOL+QROW*BCOL)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << (2.0*QROW*BROW*BCOL)/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< (4.0*QROW*BROW*BCOL+QROW*BCOL+2.0*BROW*BCOL)*DoubleSize/t1/1024/1024/1024 << endl;
  // ############################################ ----------------------------
  }
  return 0;
}


