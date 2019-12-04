#include "Ex_6lib.h"
#include <cblas.h>
#include <iostream>
using namespace std;

int main(){
	
/* 
 A: Vector x(N),Vector y(N) 
 B: Matrix A(M,N), Vector x(N)
 C: Matrix A(M,L), Matrix B(L,N)
*/ 
	int NLOOPSA, NLOOPSB, NLOOPSC;
	unsigned int NA, NB, MB, NC, MC, LC;
	NLOOPSA = 200; NA = 40000000; // exercise A
	NLOOPSB = 150; NB = 2000; MB = 40000; // exercise B
	NLOOPSC = 5; NC = 1000; MC = 1000; LC = 1200; // exercise C

		
	exA(NLOOPSA, NA);
	exB(NLOOPSB, MB, NB);
	exC(NLOOPSC, MC, NC, LC);	
				
	return 0;
}
