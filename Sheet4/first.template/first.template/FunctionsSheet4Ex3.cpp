#include "greetings.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <mpi.h>                                  // MPI
#include <string>
#include <vector>
using namespace std;


void DebugVector(vector<int> const &xin, MPI_Comm const &icomm )
{
	int myrank, numprocs; 
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
    int ProcessCalled=-1;
    
    if (0==myrank) {
        while((ProcessCalled < 0) || (ProcessCalled >= numprocs))
		{
			cout << "Please enter a valid processor number: "<< endl;
			cin >> ProcessCalled;
		}		
    }
    //Broadcast process number
    MPI_Bcast(&ProcessCalled, 1, MPI_INT, 0, icomm);
	// Print vector from called process
	if(myrank == ProcessCalled){
		cout << "Vector from process: " << myrank << " -- "; 
		for(int i=0; i < (int) xin.size(); ++i){	
			cout << xin[i] << " ";}
		cout << endl;
		cout << "-----------------------------------" << endl;
	}
	fflush(stdout);
	MPI_Barrier(icomm);
	
	// Print vectors in order
	//for(int i=0; i<numprocs; ++i ){
		//MPI_Barrier(icomm);
		//fflush(stdout);
		//if(myrank==i){
			//cout << "Vector from process: " << myrank << " -- "; 
			//for(int j=0; j < (int) xin.size(); ++j){	
			//cout << xin[j] << " ";	}
			//cout << endl;
		//}
		
	//}
	//MPI_Barrier(icomm);
    return;
    
}
