#include <cblas.h>
#include <iostream>
#include <ctime>
#include <vector>	
using namespace std;

void exA(int const NLOOPS, int const N)
{
//###########################################
//###########   InnerProduct   ################


		// Data initialization
		double re;
		int incx, incy;
    
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
	cout.precision(2);
	cout << "#######  Exercise A performance  #######"<<endl;
	cout.precision(2);
	cout << "GByte Memory:" << 2.0*N*sizeof(double)/1024/1024/1024 <<endl;
	cout.precision(6);
	cout << "sec         :" << t1 << endl;
	cout << "sec./loop   :" << t2 << endl;
	cout << "GFLOPS      :" << 2.0*N/t1/1024/1024/1024 << endl;
	cout << "GiB/s    :" << 2.0*N*sizeof(double)/t1/1024/1024/1024 << endl;
}
