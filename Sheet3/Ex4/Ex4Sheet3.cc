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

// g++ -Ofast -fopenmp Ex4Sheet3.cc && ./a.out

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

int main()
{
  cout << "Title: " << endl;
  vector<int> ThreadsUsed = {1,2,4,8,16,32,64}; //In ascending order!! The efficiency is compared wrt the efficiency of the first number of threads!!!
  const int NumberOfEvaluations = (int) ThreadsUsed.size();
  //const int InitialNumberOfThreads = 13;
  
  double EfficiencyBenchmarkMatrixVector;
  double EfficiencyBenchmarkMatrixMatrix;
  double EfficiencyBenchmarkPolynomialEvaluation;
  
  // Vector of Times (total running time)
  vector<double> TimesMatrixVector(NumberOfEvaluations);
  vector<double> TimesMatrixMatrix(NumberOfEvaluations);
  vector<double> TimesPolynomialEvaluation(NumberOfEvaluations); 
  
  // Vector of Times (total running time)
  vector<double> EfficiencyMatrixVector(NumberOfEvaluations);
  vector<double> EfficiencyMatrixMatrix(NumberOfEvaluations);
  vector<double> EfficiencyPolynomialEvaluation(NumberOfEvaluations);
  
  // Vectors of GFlops
  vector<double> GFlopsMatrixVector(NumberOfEvaluations);
  vector<double> GFlopsMatrixMatrix(NumberOfEvaluations);
  vector<double> GFlopsPolynomialEvaluation(NumberOfEvaluations);

  // Vectors of Gib/sec
  vector<double> GibPerSecMatrixVector(NumberOfEvaluations);
  vector<double> GibPerSecMatrixMatrix(NumberOfEvaluations);
  vector<double> GibPerSecPolynomialEvaluation(NumberOfEvaluations);
  
  //#################################################################### Benchmarking begins
  for(int indicator= 0; indicator<=NumberOfEvaluations; ++indicator){                 // for loop for no. of threads 
  omp_set_num_threads(ThreadsUsed[indicator]);
  double const DoubleSize = sizeof(double); // Size of double
  double tstart, t1;  // timer

  { /* */// Matrix-Vector Multiplication
  double const NLOOPSMatrixVector = 150.0; // Matrix-Vector

  if(indicator == 0){cout << "------" << endl;}
  // Test data for Matrix-Vector multiplication.
  int const MROW= 60000, NCOL= 9000;           // initialize constants

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
  TimesMatrixVector[indicator] = t1;
  if(indicator == 0){  // Setting up efficiency vector
	  EfficiencyBenchmarkMatrixVector = t1;
	  EfficiencyMatrixVector[indicator] = 1.0;}
	  else{EfficiencyMatrixVector[indicator] = EfficiencyBenchmarkMatrixVector/(ThreadsUsed[indicator]*t1);}
  //cout << "Value to use variable z: " << z[rand() % z.size()]<< endl;
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  if(indicator == 0){ cout << "Matrix-Vector Multiplication "<< endl;}
  //cout << "Total time in sec: " << t1 << endl;
  //t1 /= NLOOPSMatrixVector; // divide by number of function calls
  cout.precision(2);
  if(indicator == 0){cout << "M times N matrix: " << "M = " << MROW << ", N = " << NCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixVector << endl;}
  //cout << "Number of Threads" << "..." << endl;
  //cout << "Timing in sec. per loop: : " << t1 << endl;
  //cout << "Memory allocated (in GB): " << (MROW*NCOL + MROW + NCOL)*DoubleSize/1024/1024/1024 << endl; 
  //cout << "GFLOPS: " << 2.0*MROW*NCOL/t1/1024/1024/1024 << endl;
  GFlopsMatrixVector[indicator] = 2.0*MROW*NCOL/t1/1024/1024/1024;
  //cout << "Gib/sec: "<< (4.0*MROW*NCOL+MROW)*DoubleSize/t1/1024/1024/1024 << endl;
  GibPerSecMatrixVector[indicator] = (4.0*MROW*NCOL+MROW)*DoubleSize/t1/1024/1024/1024;
  // ############################################ ----------------------------
  }
  // ###########################################################################################
  { /* */ // Matrix-Matrix Multiplication
  if(indicator == 0){cout << "------" << endl;}
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
  // inizialization of products as zero matrices.
  for(int i; i<QROW; ++i){
	for(int j; j<BCOL; ++j){
		product[i][j]=0.0;
		product1[i][j]=0.0;
  }}
  // CheckUp of MatrixMultiplicationParallel, assert jumps if equality doesnt hold!
  if(indicator == 0)
  {
  MatrixMultiplication(Q,B,product);
  MatrixMultiplicationParallel(Q,B,product1);
  CheckEquality(product,product1);               // Verification of result
  }

  tstart = omp_get_wtime(); // start timer
  // Do Calculation
  for (int i = 0; i < NLOOPSMatrixMatrix; ++i){
    //MatrixMultiplication(Q,B,product);
    MatrixMultiplicationParallel(Q,B,product);
  }
  t1 = omp_get_wtime() - tstart;
  // End Calculation
  // ############################################ Print timing
  
  TimesMatrixMatrix[indicator] = t1;
  if(indicator == 0){  // Setting up efficiency vector
	EfficiencyBenchmarkMatrixMatrix = t1;
	EfficiencyMatrixMatrix[indicator] = 1.0;}
	else{EfficiencyMatrixMatrix[indicator] = EfficiencyBenchmarkMatrixMatrix/(ThreadsUsed[indicator]*t1);}
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  if(indicator == 0){ cout << "Matrix-Matrix Multiplication "<< endl; }
  //cout << "Total time in sec: " << t1 << endl;
  //t1 /= NLOOPSMatrixMatrix; // divide by number of function calls
  cout.precision(2);
  if(indicator == 0){ cout << "Matrix Dimensions (M times L and L times N): M = " << QROW << ", L = " << BROW << ", N = " << BCOL << endl;
  cout << "Number of Loops: " << NLOOPSMatrixMatrix << endl;}
  //cout << "Timing in sec. per loop: " << t1 << endl;
  //cout << "Memory allocated (in GB): " << (QROW*BROW + BROW*BCOL+QROW*BCOL)*DoubleSize/1024/1024/1024 << endl;
  //cout << "GFLOPS: " << (2.0*QROW*BROW*BCOL)/t1/1024/1024/1024 << endl;
  GFlopsMatrixMatrix[indicator] = (2.0*QROW*BROW*BCOL)/t1/1024/1024/1024;
  //cout << "Gib/sec: "<< (4.0*QROW*BROW*BCOL+QROW*BCOL)*DoubleSize/t1/1024/1024/1024 << endl;
  GibPerSecMatrixMatrix[indicator] = (4.0*QROW*BROW*BCOL+QROW*BCOL)*DoubleSize/t1/1024/1024/1024;
  // ############################################ ----------------------------
  }
  // ###########################################################################################
  { /* */// Polynomial Evaluation
  if(indicator == 0){cout << "------" << endl;}
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
  TimesPolynomialEvaluation[indicator] = t1;
  if(indicator == 0){  // Setting up efficiency vector
	EfficiencyBenchmarkPolynomialEvaluation = t1;
	EfficiencyPolynomialEvaluation[indicator] = 1.0;}
	else{EfficiencyPolynomialEvaluation[indicator] = EfficiencyBenchmarkPolynomialEvaluation/(ThreadsUsed[indicator]*t1);}
  //t1 /= CLOCKS_PER_SEC;     // now, t1 in seconds
  if(indicator == 0){cout << "Polynomial Evaluation "<< endl;}
  //cout << "Total time in sec: " << t1 << endl;
  //t1 /= NLOOPSPolynomialEvaluation; // divide by number of function calls
  cout.precision(2);
  if(indicator == 0){cout << "Size of vector x: " << evaluations << endl;
  cout << "Polynomial size (N) : " << polynomialSize << endl;
  cout << "Number of Loops: " << NLOOPSPolynomialEvaluation << endl;}
  //cout << "Timing in sec. per loop: " << t1 << endl;
  //cout << "Memory allocated (in GB): " << (2.0*evaluations+polynomialSize)*DoubleSize/1024/1024/1024 << endl;
  //cout << "GFLOPS: " << (3.0*polynomialSize*evaluations)/t1/1024/1024/1024 << endl;
  GFlopsPolynomialEvaluation[indicator] = (3.0*polynomialSize*evaluations)/t1/1024/1024/1024;
  //cout << "Gib/sec: "<< 7.0*evaluations*polynomialSize*DoubleSize/t1/1024/1024/1024 << endl;
  GibPerSecPolynomialEvaluation[indicator] = 7.0*evaluations*polynomialSize*DoubleSize/t1/1024/1024/1024;
  // ############################################ ----------------------------
  }
  if(indicator == 0){cout << "------" << endl;}
  }
  
  // -------------------------------  Printing Total Times in Matlab format
  cout << "TimesMatrixVector = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << TimesMatrixVector[z] << " ";}
  cout << "];" << endl;

  cout << "TimesMatrixMatrix = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << TimesMatrixMatrix[z] << " ";}
  cout << "];" << endl;
  
  cout << "TimesPolynomialEvaluation = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << TimesPolynomialEvaluation[z] << " ";}
  cout << "];" << endl; 
  cout << "------" << endl;
  
  // -------------------------------  Printing Efficiency in Matlab format
  cout << "EfficiencyMatrixVector = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << EfficiencyMatrixVector[z] << " ";}
  cout << "];" << endl;

  cout << "EfficiencyMatrixMatrix = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << EfficiencyMatrixMatrix[z] << " ";}
  cout << "];" << endl;
  
  cout << "EfficiencyPolynomialEvaluation = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << EfficiencyPolynomialEvaluation[z] << " ";}
  cout << "];" << endl; 
  cout << "------" << endl;
  
  // -------------------------------  Printing GFlops Vectors in Matlab format
  
    cout << "GFlopsMatrixVector = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GFlopsMatrixVector[z] << " ";}
  cout << "];" << endl;

  cout << "GFlopsMatrixMatrix = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GFlopsMatrixMatrix[z] << " ";}
  cout << "];" << endl;
  
  cout << "GFlopsPolynomialEvaluation = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GFlopsPolynomialEvaluation[z] << " ";}
  cout << "];" << endl;
  cout << "------" << endl;
  
  // ------------------------------- // Printing GibPerSecond Vector in Matlab format
  
  cout << "GipPerSecMatrixVector = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GibPerSecMatrixVector[z] << " ";}
  cout << "];" << endl;

  cout << "GipPerSecMatrixMatrix = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GibPerSecMatrixMatrix[z] << " ";}
  cout << "";
  cout << "];" << endl;
  
  cout << "GipPerSecPolynomialEvaluation = [ ";
  for(int z=0; z < NumberOfEvaluations; ++z)
  {cout << GibPerSecPolynomialEvaluation[z] << " ";}
  cout << "];" << endl;
  cout << "------" << endl;
  
  // ------------------------------- // Matlab plots efficiency
  cout << "NumberofThreads = [ ";
  for(int k=0; k<NumberOfEvaluations; ++k){cout << ThreadsUsed[k] << " ";}
  cout << "];" << endl;
  cout << "figure" << endl;
  cout << "subplot(3,1,1)" << endl;
  cout << "plot(NumberofThreads,EfficiencyMatrixVector)" << endl;
  cout << "subplot(3,1,2)" << endl;
  cout << "plot(NumberofThreads,EfficiencyMatrixMatrix)" << endl;
  cout << "subplot(3,1,3)" << endl;
  cout << "plot(NumberofThreads,EfficiencyPolynomialEvaluation)" << endl;
  return 0;
}
