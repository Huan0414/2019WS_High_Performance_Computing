#include <mylib.h>
using namespace std;

int main(){
	
/* 
 A: Vector x(N),Vector y(N) 
 B: Matrix A(M,N), Vector x(N)
 C: Matrix A(M,L), Matrix B(L,N)
*/ 
	int const NLOOPSA, NLOOPSB, NLOOPSC;
	unsigned int NA, N, M, L;
	NLOOPSA = 5000; NA = 1000000; // exercise A
	NLOOPSB = 50; NB = 10000; MB = 10000; // exercise B
	NLOOPSC = 10; NC = 1500; MC = 1000; LC = 1000; // exercise C

	void exA(NLOOPSA, NA);
	void exB(NLOOPSB, MB, NB);
	void exC(NLOOPSC, MC, NC, LC);	
				
	return 0;
}
