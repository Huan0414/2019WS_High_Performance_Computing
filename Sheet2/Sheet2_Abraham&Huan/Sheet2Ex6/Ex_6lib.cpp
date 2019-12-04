#include <cblas.h>
#include <iostream>
#include <ctime>
#include <vector>	
using namespace std;


//#############################################
//######        InnerProduct        ###########
void exA(int const NLOOPS, int const N)
{
		// Data initialization
		double re;
		int incx, incy;
		incx = 1;
        incy = 1;
	
	
		// This is for C instead of C++
    	// x = (double *)malloc(sizeof(double)*N);
		// y = (double *)malloc(sizeof(double)*N);
		vector<double> x(N);
		vector<double> y(N);
		for (int i=0;i<N;i++) {x[i] = (i)%219 + 1; y[i] = 1.0/x[i];}
	
		// Start calculation
		double tstart, t1, t2;
		double sk(0.0);
	
		tstart = clock(); // time started

		for(int i=0;i<NLOOPS;++i)
		{	
			re = cblas_ddot(N,x.data(),incx,y.data(),incy);
		}

		t1 = clock() - tstart;
		t1 /= CLOCKS_PER_SEC; // total time
		t2 = t1/NLOOPS; // seconds per LOOP
	
		// Check the result
		if (static_cast<unsigned int>(re) != N)
		{
			cout << "WRONG result!!\n"<<endl;
		}
	
	// Performance Evaluation
	cout << "#######  Exercise A performance  #######"<<endl;
	cout.precision(6);
	cout << "GByte Memory:" << 2.0*N*sizeof(double)/1024/1024/1024 <<endl;
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2.0*N/t2/1024/1024/1024 << endl;
	cout << "GiB/s       :" << 4.0*N*sizeof(double)/t2/1024/1024/1024 << endl;
}

//#############################################
//#############   Matrix Vector   #############	
void exB(int const NLOOPS,int const M, int const N)
{	
	// Data Initialization	
	CBLAS_LAYOUT Layout;
	CBLAS_TRANSPOSE transa;
	
    double alpha, beta;
    int LDA,incx,incy;
    
    Layout = CblasColMajor;
    transa = CblasNoTrans;	
	LDA = M;
    incx = 1;
    incy = 1;
	alpha = 1;
	beta = 0;
	
    vector<double> A(M*N);
	vector<double> x(N);
	vector<double> y(M);

    
    for (int i=0;i<M;i++) 
    {
		for (int j=0;j<N;j++)
			{
				A[i+j*M]= (i+j)%219 + 1;  
				if(i==17) {x[j] = 1.0/A[17+j*M];y[j] = 0;}
			}
	}
	
	// Do Calculation	
	double tstart, t1, t2;
	double sk(0.0);
	tstart = clock(); // time started
	
	for(int i=0;i<NLOOPS;++i)
	{	
		cblas_dgemv(Layout,transa,M,N,alpha,A.data(),LDA,x.data(),incx,beta,y.data(),incy);
	}
	
	t1 = clock() - tstart;
	t1 /= CLOCKS_PER_SEC; // total time
	t2 = t1/NLOOPS; // seconds per LOOP
	
	// Performance Evaluation
	cout << "#######  Exercise B performance  #######"<<endl;
	cout.precision(6);
	cout << "GByte Memory:" << 1.0*(N*M + N + M)*sizeof(double)/1024/1024/1024 <<endl;
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2.0*N*M/t2/1024/1024/1024 << endl;
	cout << "GiB/s       :" << (4.0*M*N + M)*sizeof(double)/t2/1024/1024/1024 << endl;
}


//#############################################
//#############   Matrix Matrix   #############	
void exC(int const NLOOPS, int const M,int const N,int const K)
{
	// Data initialization
	CBLAS_LAYOUT Layout;
	CBLAS_TRANSPOSE transa;
	CBLAS_TRANSPOSE transb;
	
    double alpha, beta;
    int LDA,LDB,LDC;

    Layout = CblasColMajor;
    transa = CblasNoTrans;
    transb = CblasNoTrans;

    LDA = M;
    LDB = K;
    LDC = M;
	alpha = 1;
	beta = 0;

    double *A = new double[M * K]; // define as static array
    double *B = new double[K * N];
    double *C = new double[M * N];

    for (int i=0;i<M;i++) {
		for(int j=0;j<K;j++) {A[i+j*M]= rand()%10 + 1;}}  // Matrix A
	for(int j=0;j<K;j++) {
		for (int k=0;k<N;k++) {B[j+k*K] = rand()%10 + 1;}} // Matrix B	
	for (int i=0;i<M;i++) {
		for (int k=0;k<N;k++) {C[i+k*M] = 0;}} // Matrix C is filled with zero

    //  Do Calculation
    double tstart, t1, t2;
	double sk(0.0);
    tstart = clock(); // time started
	
	for(int i=0;i<NLOOPS;++i)
	{	
    cblas_dgemm(Layout,transa,transb,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC);
	}
	
	t1 = clock() - tstart;
	t1 /= CLOCKS_PER_SEC; // total time
	t2 = t1/NLOOPS; // seconds per LOOP
	
	// Performance Evaluation
	cout << "#######  Exercise C performance  #######"<<endl;
	cout.precision(6);
	cout << "GByte Memory:" << 1.0*(M*K + K*N + M*N)*sizeof(double)/1024/1024/1024 <<endl;
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2.0*N*M*K/t2/1024/1024/1024 << endl;
	cout << "GiB/s       :" << (4.0*M*K*N + M*N)*sizeof(double)/t2/1024/1024/1024 << endl;
	
	delete [] A;
	delete [] B;
	delete [] C;

}

