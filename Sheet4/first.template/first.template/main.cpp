//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "greetings.h"
#include <iostream>                              // MPI
#include <mpi.h>
#include <vector>
#include <numeric> 								// iota
#include <cassert>
using namespace std;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);    
    MPI_Comm icomm = MPI_COMM_WORLD;
	//int MyVectorSize = 4;
    int myrank, numprocs;
    
    //numprocs = 1;    // delete this line when uncommenting the next line
    MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
	//vector<int> vec(MyVectorSize); 
	//iota(vec.begin(),vec.end(),myrank);						// Vector Initialization
		
	//cout << "Greetings from process " << myrank << endl;
	//MPI_Barrier(icomm); 
	
    if (0==myrank) {
        cout << "\n There are " << numprocs << " processes running.\n \n";	
    }
    
    //DebugVector(vec,icomm);
	//PrintVectorsInOrder(vec,icomm);
    //greetings(icomm);
    //greetings_cpp(icomm);
    //fflush(stdout);
	//MPI_Barrier(icomm);
	
	//// ########################## Dot Product #############################
	//{if(myrank == 0) cout << "##################################################### Ex 7 -- Dot Product" << endl ;
	//double innerProduct = 0.0;
	//size_t const N = 10;								  // VectorSize DocProduct
	//vector<double> vec1(N,0);
	//vector<double> vec2(N,0);
	
	//// Vectors Initialization
	//if(myrank == 0)
	//{
		////Data initialization
		//for (unsigned int i = 0; i < N; ++i)
		//{
			//vec1[i] = i + 1.0;
			////cout << x[i] << " ";
			//vec2[i] = 1.0 / vec1[i];
		//}
	//}
	
	//par_scalar(innerProduct,vec1,vec2,icomm);				// Calculation of ScalarProduct in MPI
	
	//if(myrank == 0) cout << "Process: " << myrank << " -- Inner Product: " << innerProduct << endl; 
	////cout << "Process: " << myrank << " -- Inner Product: " << innerProduct << endl;
	//}

	////############################# E7 Task3 MinMaxExchange ###########################
	if(myrank == 0) cout << "##################################################### Ex 7 MinMaxExchange" << endl ;
			
	//vector<double> x;	
	//if (0==myrank) x = {1.0,2.0,5.0,-59.0,-102.0,-23.0,-50.0, 90.0, 30.0};
	// Scatter vectors to all the processors
	//int N = x.size();
	
	//#################################################################################
	//if(myrank == 0) cout << "Scatter Count value: " << ScatterCount << endl;
	//vector<double> localVec(ScatterCount);	
	//MPI_Scatter(&x[0],ScatterCount,MPI_DOUBLE,&localVec[0],ScatterCount,MPI_DOUBLE,0,icomm);	// Why scatter is failing?
	//MPI_Barrier(icomm);	
	
	int VectorSize = 4;
	vector<double> localVec(VectorSize);
	
	if (myrank == 0) localVec = {4.0,2.0,40.0,1000.0};
	if (myrank == 1) localVec = {5.0,-59.0,-1000.0,-1000.0};
	if (myrank == 2) localVec = {1000.0,-23.0,-1000.0,-78.0};
	if (myrank == 3) localVec = {-50.0, 1000.0,21.0,0.0};

	
	// Find Global Min Max
	MinMaxWithoutReduc(localVec,icomm);
	
	
	////############################# E8 Task3 #########################################
	//if(myrank == 0) cout << "##################################################### Ex 8" << endl ;
	//MPI_Barrier(icomm);
	//int SizeVector = 20;
	//vector<double> x(SizeVector);
	//vector<double> exchanged(SizeVector);
	////vector<double> transference(5);
	
	////// Initialization of x and transference vectors
	//for(int i=0; i<SizeVector; ++i) x[i] = myrank*100.0 + (i%5)*10 + i;
	////for(int j=0; j<5; ++j) transference[j] = myrank*5.0 + j;

	////// Print x and transference vectors
	//if(myrank == 0) cout << "BEFORE transference" << endl;
	//PrintVectorsInOrder(x,icomm);
	////PrintVectorsInOrder(transference,icomm);
	
	//MPI_Alltoall(&x[0], 5, MPI_DOUBLE, &exchanged[0], 5, MPI_DOUBLE, icomm);
	////MPI_Alltoall(MPI_IN_PLACE, 5, MPI_DOUBLE, &x[0], 5, MPI_DOUBLE, icomm);
	
	//if(myrank == 0) cout << "AFTER transference" << endl;
	//PrintVectorsInOrder(exchanged,icomm);
		
    MPI_Finalize();

    return 0;
}


