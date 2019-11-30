#include <cstdio>
#include <cstdlib>
#include <cblas.h>
#include <iostream>
using namespace std;

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
	cout.precision(2);
	cout << "#######  Exercise C performance  #######"<<endl;
	cout.precision(2);
	cout << "GByte Memory:" << (M*K + K*N + M*N)*sizeof(double)/1024/1024/1024 <<endl;
	cout.precision(6);
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2*N*M*K/t1/1024/1024/1024 << endl;
	cout << "GiB/s       :" << (M*K + K*N + M*N)*sizeof(double)/t1/1024/1024/1024 << endl;
	
	delete [] A;
	delete [] B;
	delete [] C;

}
