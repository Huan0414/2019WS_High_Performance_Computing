#include <ctime>
#include <iostream>
#include <omp.h>            // OpenMP
#include <assert.h>
#include <cmath>
#include <vector>
#include "mylib.h"
using namespace std;

int main()
{
    double const NLOOPSPolynomialEvaluation = 10.0;
	int polynomialSize = 200000;
	int evaluations = 5000;
	vector<double> coefficients(polynomialSize);
	vector<double> x(evaluations);
	
	// test data
	for(int i=0; i<evaluations; i++){
		//x[i] = 1;
		x[i] = i%219+1;
	}

	for(int i=0; i<polynomialSize; i++){
		//coefficients[i] = i+1;
		coefficients[i] = (i)%219 + 1;
	}
	
	// Do Calculation
	vector<double> result(evaluations, 0.0);
	double const DoubleSize = sizeof(double); // Size of double
	double tstart, t1;  // timer
	
	// Clock starts and Calculation begins
	tstart = omp_get_wtime();
	for (int i=0; i < NLOOPSPolynomialEvaluation; ++i){
		//PolynomialEvaluation(coefficients,x,result);  // without openmp
		//PolynomialEvaluation1(coefficients,x,result);  // parallize outer loop
		PolynomialEvaluation2(coefficients,x,result); 	// parallize inner loop
	}
	// End Calculation
	// ############################################ Print timing
	t1 = omp_get_wtime() - tstart;
  
	cout << "Total time in sec: " << t1 << endl;
	t1 /= NLOOPSPolynomialEvaluation; // divide by number of function calls
    cout << "Timing in sec. per loop: " << t1 << endl;
	
	//for(int j = 0;j < evaluations;j++)
	//{cout << result[j]<<" ";}
	
	
	return 0;
}

//###################################################################
//polynomial function, parallize outer loop
void PolynomialEvaluation(vector<double> const &a, vector<double> const &x, vector<double> &result){
    int Xsize = x.size();
    int Asize = a.size();
    int Rsize = result.size();
    assert(Rsize == Xsize);
    
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

//###################################################################
//polynomial function, parallize outer loop
void PolynomialEvaluation1(vector<double> const &a, vector<double> const &x, vector<double> &result){
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

//###################################################################
//polynomial function, parallize inner loop
void PolynomialEvaluation2(vector<double> const &a, vector<double> const &x, vector<double> &result){
    int Xsize = x.size();
    int Asize = a.size();
    int Rsize = result.size();
    assert(Rsize == Xsize);
	    
            
    for (int j = 0; j< Xsize; j++)
    {
        double sum = 0.0;
        
        #pragma omp parallel default(none) shared(a,x,Xsize,Asize,j) reduction(+:sum)
		{
			
			int nthrd = omp_get_num_threads();
			int tid = omp_get_thread_num();
			int portionsize = 1.0*Asize/nthrd; // divide calculation to all threads
			int sizestart = tid*portionsize ; // start number of each thread
			int sizeend = min(sizestart + portionsize , Asize + 1); // end number of each thread
			
			double x_iterative = pow(x[j],sizestart); //define initialized x^k for each thread
			for (int i = sizestart; i < sizeend; i++)
			{
				sum += a[i]*x_iterative;
				x_iterative *= x[j]; 
			}
		}
		result[j] = sum;
	}
	
}
