#include <lapacke.h>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;

int main(){
	
	double *A, *B;
	int *IPIV;
	int Layout, n, nrhs, LDA,LDB, INFO; 
	char trans = 'N'; //No transpose;
	Layout = LAPACK_COL_MAJOR;
    INFO = 3;
	
	double tstart, t1, t2;
	double sk(0.0);

	for(int countn=0;countn<3;countn++) // n = {500, 1000, 2000}
	{  
		cout << endl;
	   if(countn==0) {n = 500; cout << "########  When n is 500  ######"<< endl;}
	   if(countn==1) {n = 1000; cout << "########  When n is 1000  ######"<< endl;}
	   if(countn==2) {n = 2000; cout << "########  When n is 2000  ######"<< endl;}
	    LDA = n;
        LDB = n;
        nrhs = 0;
		
		for (int countnr = 0;countnr < 9 ;countnr++)	
		{			
			nrhs += 50;
			
		// Initialize matrix A(n*n) and right hand side matrix B(n*nrhs)
			double *A = new double[n*n];
			double *B = new double[n*nrhs];
			int *IPIV = new int[n];
		/*
		* ============================================
		* 1. For test, n=3, nrhs = 1, output: 11 -9  1
		* ============================================
		*/
		//A[0]=1; A[1]=2; A[2]=3;A[3]=2;A[4]=3;A[5]=4;A[6]=3;A[7]=4;A[8]=1; 
		//B[0]=-4; B[1]=-1;B[2]=-2;	
		
		/*
		* ==========================================
		* 2. Set Matrix A and B randomly
		* ==========================================
		*/		
			for (int i=0;i<n;i++) 
			{
				for (int j=0;j<n;j++)
					{A[i+j*n]= rand()%10 + 1;}
			}
			for (int i=0;i<n;i++) 
			{
				for (int j=0;j<nrhs;j++)
					{B[i+j*n]= rand()%10 + 1;}
			}

			// Do calculate
			tstart = clock(); // time started	
			LAPACK_dgetrf(&n,&n,A,&LDA,IPIV,&INFO);

			if(INFO) {cout << "An error occured" <<endl;}
			else {
			//		cout << "Got the LU factorization, starting solving equation..."<<endl;
					LAPACK_dgetrs(&trans,&n,&nrhs,A,&LDA,IPIV,B,&LDB,&INFO);
						if(INFO) {cout<<"An error occuring in dgestrs"<<endl;}}
			//		else{
			//			cout << "The solution is:" ;
			//			for (int i=0;i<n;i++) {cout << B[i]<<" ";}}}
			
			t1 = clock() - tstart;
			t1 /= CLOCKS_PER_SEC; // total time
			t2 = t1/nrhs; // seconds per right hand side
			cout << "sec./nrhs = "<< t2 << " nrhs = "<<nrhs<<endl;
		}
	}
	
	delete [] A;
	delete [] B;
	delete [] IPIV;
				
	return 0;
}

