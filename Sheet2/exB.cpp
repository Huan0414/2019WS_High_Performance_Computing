#include <cblas.h>
#include <iostream>
#include <ctime>
#include <vector>
using namespace std;	

void exB(int const NLOOPS,int const M, int const N)
{	
//##########################################
//#############   Matrix Vector   #############	
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
	cout.precision(2);
	cout << "#######  Exercise B performance  #######"<<endl;
	cout.precision(2);
	cout << "GByte Memory:" << (N*M + N + M)*sizeof(double)/1024/1024/1024 <<endl;
	cout.precision(6);
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2*N*M/t1/1024/1024/1024 << endl;
	cout << "GiB/s       :" << (N*M + N + M)*sizeof(double)/t1/1024/1024/1024 << endl;
}
