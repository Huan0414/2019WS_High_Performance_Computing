#include <cstdlib>          // atoi()
#include <cstring>          // strncmp()
#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <sstream>
#include <string>
#include <assert.h>
#include <cmath>
#include <vector>

using namespace std;

double InnerProduct(vector<double> const &x, vector<double> const &y, double &sum)
{
    assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    //double sum = 0.0;
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
  for(int i=0; i<RowsA; ++i){
    for(int j=0; j<ColumnsA; ++j){
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
      assert(B[j][i] == A[i][j]); }}
}

void MatrixMultiplicationParallel( std::vector<std::vector<double>> const &a,  std::vector<std::vector<double>> const &b,
                          std::vector<std::vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   assert(columnA == rowB); // Checks if dimensions match
   std::vector<std::vector<double>> TransposeB(columnB,std::vector<double>(rowB));
   TransposeMatrix(b,TransposeB);
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   #pragma omp parallel for default(none) shared(a,product,rowA,columnB,TransposeB) collapse(2)
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

int main(){
  double const DoubleSize = sizeof(double); // Size of double
  double tstart, t1;  // timer

  { /* */// Matrix-Vector Multiplication
  double const NLOOPSMatrixVector = 150.0; // Matrix-Vector

  cout << "------" << endl;
  // Test data for Matrix-Vector multiplication.
  int const MROW= 6000, NCOL= 9000;           // initialize constants

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
  tstart = omp_get_wtime(); // start timer
  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixVector; ++i){
     MatrixVectorMultiplication(A,x,z);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = omp_get_wtime() - tstart;
  cout << "Value to use variable z: " << z[rand() % z.size()]<< endl;
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  cout << "Matrix-Vector Multiplication "<< endl;
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixVector; // divide by number of function calls
  cout.precision(2);
  cout << "M times N matrix: " << "M = " << MROW << ", N = " << NCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixVector << endl;
  cout << "Number of Threads" << "..." << endl;
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

  int const QROW= 2000, BROW = 2000, BCOL= 2000;
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
  vector<vector<double>> product1(QROW, vector<double>(BCOL));
  MatrixMultiplication(Q,B,product);
  MatrixMultiplicationParallel(Q,B,product1);
  CheckEquality(product,product1);               // Verification of result
  tstart = omp_get_wtime(); // start timer

  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixMatrix; ++i){
    MatrixMultiplicationParallel(Q,B,product);
    //MatrixMultiplicationParallel(Q,B,product);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = omp_get_wtime() - tstart;
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
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
  // Clock starts and Calculation begins
  tstart = omp_get_wtime();
  for (int i=0; i < NLOOPSPolynomialEvaluation; ++i){
    PolynomialEvaluation(coefficients,x,result);
  }
  // End Calculation
  // ############################################ Print timing
  t1 = omp_get_wtime() - tstart;
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
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
  return 0;
}
