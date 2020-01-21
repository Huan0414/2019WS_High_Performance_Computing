//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "greetings.h"
#include <iostream>                               // MPI
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    MPI_Comm icomm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    int myrank, numprocs;
    //numprocs = 1;    // delete this line when uncommenting the next line
    MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);

    if (0==myrank) {
      cout << "\n There are " << numprocs << " processes running.\n \n";
    }
	
	vector<double> x;	
	if (0==myrank) {
	x = {1.0,2.0,5.0,-59.0,-102.0,-23.0,-50.0};
	}

	MinMaxWithoutReduc(x,icomm);
	//MinMaxWithReduc(x,icomm);
    //greetings(icomm);
    //greetings_cpp(icomm);

    //if (0==myrank) cout << endl;

    MPI_Finalize();

    return 0;
}

void MinMaxWithoutReduc(vector<double> const &x, MPI_Comm const icomm)
{
	int myrank, numprocs;
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);              // #MPI processes
	vector<double> globalMinVec(numprocs,0.0), globalMaxVec(numprocs,0.0);


// Scatter vectors to all the processors
	int N = x.size();
	int ScatterCount = max(N/numprocs,1);
	vector<double> localVec(ScatterCount,0.0);	
	MPI_Scatter(&x[0],ScatterCount,MPI_DOUBLE,&localVec[0],ScatterCount,MPI_DOUBLE,0,icomm);	
	
// Each process find their min and max number
	double localMin =  *min_element (localVec.begin(), localVec.end()); 
	//cout <<"Processor "<< myrank << " has localMin" << localMin<<endl;
	double localMax =  *max_element (localVec.begin(), localVec.end());
	
// Processor 0 gather all the result to a vector	
	MPI_Gather(&localMin,1,MPI_DOUBLE,&globalMinVec[0],1,MPI_DOUBLE,0,icomm);
	MPI_Gather(&localMax,1,MPI_DOUBLE,&globalMaxVec[0],1,MPI_DOUBLE,0,icomm);
    
// Processor 0 deal with elements left and find the global min, max number
	double globalMin, globalMax;
	if (0==myrank) {
		int leftover = N - ScatterCount*numprocs;
		globalMin =  *min_element (globalMinVec.begin(), globalMinVec.end());
		globalMax =  *max_element (globalMaxVec.begin(), globalMaxVec.end());
		
		if(leftover>0)
		{
			double leftMin;
			double leftMax;
			
			//leftMin = *min_element (x.begin()+N-leftover, x.end());
			//leftMax = *max_element (x.begin()+N-leftover, x.end());
			leftMin = *min_element (x.begin()+N-leftover, x.end());
			leftMax = *max_element (x.begin()+N-leftover, x.end());
			
			globalMin = min(globalMin,leftMin);
			globalMax = max(globalMax,leftMax);
		}
		cout << "The max number in the vector is: " << globalMax << endl;
		cout << "The min number in the vector is: " << globalMin << endl;
		}
		
	MPI_Bcast(&globalMax,1,MPI_DOUBLE,0,icomm);
	MPI_Bcast(&globalMin,1,MPI_DOUBLE,0,icomm);
	//cout <<"In Processor "<< myrank << " the min number is: " << globalMin << endl;
	
	for(int i = 0; i < ScatterCount;i++){
		if(localVec[i]==globalMin){localVec[i]=globalMax;}
		else if(localVec[i]==globalMax){localVec[i]=globalMin;}
	}
	
	vector<double> NewX(N,0.0);
	MPI_Gather(&localVec[0],ScatterCount,MPI_DOUBLE,&NewX[0],ScatterCount,MPI_DOUBLE,0,icomm);
	
	if(0==myrank){
		int leftover = N - ScatterCount*numprocs;
		if(leftover>0) {
			for(int i = 0; i < leftover;i++){
				int newIndex = numprocs*ScatterCount+i;
				NewX[newIndex] =  x[newIndex];
				if(NewX[newIndex]==globalMin) {NewX[newIndex]=globalMax;}
				else if(NewX[newIndex]==globalMax) {NewX[newIndex]=globalMin;}
		}
		}
		cout<<"Final x is: "<<endl; for(int i = 0; i < N; i++) {cout<<" "<<NewX[i];}
		cout << endl;
	}
		
}

// Put this into head file
void MinMaxWithoutReduc(std::vector<double> const &x, MPI_Comm const icomm);
