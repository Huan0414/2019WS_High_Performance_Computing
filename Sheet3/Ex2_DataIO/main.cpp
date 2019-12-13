#include <iostream>
#include <vector>
#include <string>
#include <ctime>
using namespace std;

#include "file_io.h"
#include "mylib.h"

// g++ -O3 main.cpp mylib.cpp file_io.cpp && ./a.out

int main()
{
 //###################################################################
 // Initialization
    cout << "File einlesen." << endl;
    const string name("data_1.txt");        // name of input file
    const string name2("out_1.txt");        // name of output file
	int NLOOPS = 20000;
    
//####################################################################    
// read, calculate and write
	double tstart, t1, t2;  // timer
	tstart = clock();
	
	vector<double> a;
    vector<double> b(6);
    size_t const N = a.size();
	
	for(int i=0; i< NLOOPS; i++)
	{
    read_vector_from_file(name, a);	 // read vector from file	
    cal_min_max_mean(a, b);			 // get min,max and mean values    		
    write_vector_to_file(name2, b);	 // get min,max and mean values
	}
	
    t1 = clock() - tstart;
    t1 /= CLOCKS_PER_SEC;	// t1 in seconds
	t2 = t1/NLOOPS;			// Time for each loop
//####################################################################    
// Performance evaluation
    cout << endl;
    cout.precision(4);
    cout << "Total time	:" << t1 << endl;
    cout << "Time per loop	:" << t2 << endl;
    //cout << "GFLOPS		:" << (4.0*N + 2) / t2 / 1024 / 1024 / 1024 << endl;
    //cout << "GiByte/s	:" <<  (N + 6) / t1 / 1024 / 1024 / 1024 * sizeof(double) << endl;
    
//###################################################################    
// print out
	cout << "Final result is	: " ;
    for (unsigned int k=0; k<b.size(); ++k)
    {
        cout << "  " << b.at(k);
    }
    cout << endl;

    return 0;
}