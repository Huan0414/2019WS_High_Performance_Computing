#include<iostream>
#include<vector>
#include<assert.h>
#include<ctime>
#include<cmath>
#include <cstdlib>

using namespace std;



// Function Inner Product
void InnerProduct(vector<double> const &vecOne, vector<double> const &vecTwo, double &product)
{
	assert(vecOne.size() == vecTwo.size());
    product = 0.0;
    for (int i = 0; (unsigned)i <= vecOne.size()-1; i++){
      product = product + (vecOne[i])*(vecTwo[i]);
    }
    //return product;
}

// L2 Norm of vector
void LTwoNorm(vector<double> const &vec, double &product)
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
void MatrixVectorMultiplication(vector<vector<double>> const &A, vector<double> const &x, vector<double> &product)
{
    int const mrows = A.size();         // #matrix rows

    for (int i = 0; i < (unsigned)mrows; ++i) {
      InnerProduct(A[i],x,product[i]);                    // inner product of row i with vector x
      //cout<< product[i] << " ";
    }
}

// Function TransposeMatrix
void TransposeMatrix(vector<vector<double>> const &A, vector<vector<double>> &Transpose){
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
void MatrixMultiplication( vector<vector<double>> const &a,  vector<vector<double>> const &b,
                          vector<vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   assert(columnA == rowB); // Checks if dimensions match
   vector<vector<double>> TransposeB(columnB,vector<double>(rowB));
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

// Function MatrixMultiplication Column-wise access
void MatrixMultiplicationColumnwiseAccess( vector<vector<double>> const &a,
                                          vector<vector<double>> const &b, vector<vector<double>> &product) {
   int rowA=a.size(), columnA=a[0].size(), rowB=b.size(), columnB=b[0].size();
   //vector<vector<double>> product(rowA, vector<double>(columnB));
   assert(columnA == rowB); // Checks if dimensions match
   for(int i=0; i<rowA; ++i){
     for(int j=0; j<columnB; ++j){
       product[i][j] = 0;
       for(int k=0; k<columnA; ++k){
         product[i][j]+=a[i][k]*b[k][j];
   }}}
   //return C;
   //cout<<"The product is:"<<endl;
   //for(int i=0; i<rowA; ++i) {
   //  for(int j=0; j<columnB; ++j){
   //   cout<<product[i][j]<<" ";}
   //  cout<<endl;}
}

void PolynomialEvaluation(vector<double> const &a, vector<double> const &x, vector<double> &result){
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

// ############################################### Main
int main()
{
   // Chose a value such that the benchmark runs at least 10 sec.
  double const NLOOPSInnerProd = 200.0; // InnerProduct
  double const NLOOPSL2 = 200.0; // L2Norm
  double const NLOOPSMatrixVector = 200.0; // Matrix-Vector
  double const NLOOPSMatrixMatrix = 1.0; // Matrix-Matrix
  double const NLOOPSPolynomialEvaluation = 10.0;
  double const DoubleSize = sizeof(double);

  double tstart, t1;  // timer

  /* *///############################################ Inner Product
  // Test data for InnerProduct
  {int const VectorSize = 6000000;
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
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSInnerProd;           // divide by number of function calls
  cout << endl;
  cout.precision(2);
  cout << "InnerProduct "<< endl;
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Memory allocated: " << 2.0*VectorSize*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*VectorSize/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 2.0*VectorSize*DoubleSize/t1/1024/1024/1024 << endl;
  }

  //############################################# -------------

  /* */// ############################################ Matrix-Vector Multiplication
  cout << "------" << endl;
  // Test data for Matrix-Vector multiplication.
  {int const MROW=10000, NCOL=900;           // initialize constants

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
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixVector; // divide by number of function calls
  cout << endl;
  cout.precision(2);
  cout << "Matrix-Vector Multiplication "<< endl;
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Memory allocated (in GB): " << (MROW*NCOL + MROW + NCOL)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*MROW*NCOL/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< (MROW*NCOL+NCOL+MROW)*DoubleSize/t1/1024/1024/1024 << endl;
  }
  // ############################################ ----------------------------

  /* */ // ############################################ Matrix-Matrix Multiplication
  cout << "------" << endl;
  {int const QROW= 1000, BROW = 1000, BCOL= 1000;
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
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSMatrixMatrix; // divide by number of function calls
  cout << endl;
  cout.precision(2);
  cout << "Matrix-Matrix Multiplication "<< endl;
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Memory allocated (in GB): " << (QROW*BROW + 2.0*BROW*BCOL+QROW*BCOL)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << (2.0*QROW*BROW*BCOL)/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< (QROW*BROW+3.0*BROW*BCOL+QROW*BCOL)*DoubleSize/t1/1024/1024/1024 << endl;
  }
  // ############################################ ----------------------------

  /* */// ############################################ Polynomial Evaluation
  cout << "------" << endl;
  // test data
  {int polynomialSize = 9000000;
  int evaluations = 25;
  vector<double> coefficients(polynomialSize);
  vector<double> x(evaluations);

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
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSPolynomialEvaluation; // divide by number of function calls
  cout << endl;
  cout.precision(2);
  cout << "Polynomial Evaluation "<< endl;
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Memory allocated (in GB): " << (2.0*evaluations+polynomialSize)*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << (3.0*polynomialSize*evaluations)/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< evaluations*(evaluations+3.0*polynomialSize)*DoubleSize/t1/1024/1024/1024 << endl;
  }

  // ############################################ ----------------------------

  /* *///############################################ L2 Norm
  cout << "------" << endl;
  // Test data for L2Norm
  {int const VectorSize = 9000000;
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
  cout << "Total time in sec: " << t1 << endl;
  t1 /= NLOOPSL2;           // divide by number of function calls
  cout << endl;
  cout.precision(2);
  cout << "L2 Norm "<< endl;
  cout << "Timing in sec. : " << t1 << endl;
  cout << "Memory allocated: " << VectorSize*DoubleSize/1024/1024/1024 << endl;
  cout << "GFLOPS: " << 2.0*VectorSize/t1/1024/1024/1024 << endl;
  cout << "Gib/sec: "<< 2.0*VectorSize*DoubleSize/t1/1024/1024/1024 << endl;}

  //############################################# -------------


  return 0;
}
