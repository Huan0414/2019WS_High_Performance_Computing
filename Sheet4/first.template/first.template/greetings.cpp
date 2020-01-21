#include "greetings.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <mpi.h>                                  // MPI
#include <string>
#include <vector>
#include <algorithm>    						// std::min
using namespace std;

// see  http://www.open-mpi.org/doc/current
// for details on MPI functions

void DebugVector(vector<double> const &xin, MPI_Comm const &icomm )
{
	MPI_Barrier(icomm);
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

void PrintVectorsInOrder(vector<double> const &xin, MPI_Comm const &icomm ){
	int myrank, numprocs; 
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
	vector<double> data(xin.size());
	int ierr;
	// prcessor 0 prints vector
	if(0==myrank){
	cout << "Vector from process: " << myrank << " -- "; 
	for(int i=0; i < (int) xin.size(); ++i){	
		cout << xin[i] << " ";}
	cout << endl;
	
	// Other processors print their vectors
	for (int i = 1; i < numprocs; ++i) 
	{
		MPI_Status stat;
        stat.MPI_ERROR = 0;            // M U S T   be initialized!!
		ierr = MPI_Recv(&data[0], xin.size(), MPI_DOUBLE, i, MPI_ANY_TAG, icomm, &stat); // MPI_ANY_SOURCE to i to order
		assert(0==ierr);
		
		cout << "Vector from process: " << i << " -- "; 
		for(int j=0; j < (int) xin.size(); ++j){	
			cout << data[j] << " ";}
		cout << endl;				
	   }
	}
	else
	{
		int dest = 0;
        ierr = MPI_Send(&xin[0], xin.size(), MPI_DOUBLE, dest, myrank, icomm);
        assert(0==ierr);
	}	
	return;
	}

void scalarProduct(double &sum, vector<double> const &x, vector<double> const &y, MPI_Comm const &icomm)
{
	int myrank, numprocs;
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);              // #MPI processes
	
	assert(x.size() == y.size()); // switch off via compile flag: -DNDEBUG
    size_t const N = x.size();
    sum = 0.0;
    #pragma omp parallel for default(none) shared(x,y,N) reduction(+:sum) 
    for (size_t i = 0; i < N; ++i)
    {
        sum += x[i] * y[i];
        //sum += exp(x[i])*log(y[i]);
    }
    return;
}

void par_scalar(double &innerProduct, vector<double> const &vec1, vector<double> const &vec2, MPI_Comm const &icomm)
{
	int myrank, numprocs;
	MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
	size_t const N = vec1.size();			// VectorSize DocProduct
	assert(N == vec2.size());
	double localInnerProduct = 0.0;								  
	size_t const sentElements = max( (int) N/numprocs,0);
	vector<double> localVec1(sentElements), localVec2(sentElements);	// Local vectors
	
	// Scatter of data among processors
	MPI_Scatter(&vec1[0], sentElements, MPI_DOUBLE, &localVec1[0], \
	sentElements, MPI_DOUBLE, 0, icomm);
	MPI_Scatter(&vec2[0], sentElements, MPI_DOUBLE, &localVec2[0], \
	sentElements, MPI_DOUBLE, 0, icomm);
	//fflush(stdout);
	MPI_Barrier(icomm);	
	
	// Calculation of local inner products
	if(myrank==0){
		//cout << "Sent elements: " << sentElements << endl ;
		scalarProduct(localInnerProduct,localVec1,localVec2, icomm);   // Innner Product in rank 0		
		// Lefovers Handling
		int leftOver = N - sentElements*numprocs;          // elements remaining?
		double innerProductLeftOver = 0.0;
		if(leftOver>0){		
			vector<double> LeftOverVector1(leftOver,0), LeftOverVector2(leftOver,0);
			for(int j=0; j< (int)LeftOverVector1.size(); ++j){
				LeftOverVector1[j] = vec1[sentElements*numprocs + j];
				LeftOverVector2[j] = vec2[sentElements*numprocs + j];
			}
			scalarProduct(innerProductLeftOver, LeftOverVector1, LeftOverVector2, icomm);
			localInnerProduct += innerProductLeftOver;
		}
		
		//cout << "Processor: " << myrank << "-- LocalInnerProduct: " << localInnerProduct << endl;
		//fflush(stdout);	
	}
	else{
		scalarProduct(localInnerProduct,localVec1,localVec2, icomm);  // scalar product of vectors outside rank 0
		//cout << "Processor: " << myrank << " -- LocalInnerProduct: " << localInnerProduct << endl;
		//fflush(stdout);
	}
	
	// Reduction of inner products
	MPI_Barrier(icomm);
	MPI_Allreduce(&localInnerProduct,&innerProduct,1,MPI_DOUBLE, MPI_SUM, icomm);
	//MPI_Reduce(&localInnerProduct,&innerProduct,1,MPI_DOUBLE,MPI_SUM, 0, icomm);
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

void greetings(MPI_Comm const &icomm)
{
    int myrank, numprocs;
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);              // #MPI processes
    char *name  = new char [MPI_MAX_PROCESSOR_NAME],
         *chbuf = new char [MPI_MAX_PROCESSOR_NAME];

    int reslen, ierr;
    MPI_Get_processor_name( name, &reslen);

    if (0==myrank) {
        cout << "   " << myrank << " runs on " << name << endl;
        for (int i = 1; i < numprocs; ++i) {
            MPI_Status stat;
            stat.MPI_ERROR = 0;            // M U S T   be initialized!!

		    ierr = MPI_Recv(chbuf, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, \
							i, MPI_ANY_TAG, icomm, &stat); // MPI_ANY_SOURCE to i to order
            assert(0==ierr);

            cout << "   " << stat.MPI_SOURCE << " runs on " << chbuf;
            int count;
            MPI_Get_count(&stat, MPI_CHAR, &count); // size of received data
            cout << " (length: " << count << " )" << endl;
            //  stat.Get_error()	// Error code
        }
    }
    else {
        int dest = 0;
        ierr = MPI_Send(name, strlen(name) + 1, MPI_CHAR, dest, myrank, icomm);
        assert(0==ierr);
    }
    delete [] chbuf;
    delete [] name;
    return;
}


void greetings_cpp(MPI_Comm const &icomm)
{
    int myrank, numprocs;
    MPI_Comm_rank(icomm, &myrank);                // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);              // #MPI processes
    string name(MPI_MAX_PROCESSOR_NAME,'#'),      // C++
           recvbuf(MPI_MAX_PROCESSOR_NAME,'#');   // C++: receive buffer, don't change size

    int reslen, ierr;
    MPI_Get_processor_name(name.data(), &reslen);
    name.resize(reslen);                          // C++

    if (0==myrank) {
        cout << "   " << myrank << " runs on " << name << endl;
        for (int i = 1; i < numprocs; ++i) {
            MPI_Status stat;
            stat.MPI_ERROR = 0;            // M U S T   be initialized!!

            ierr = MPI_Recv(recvbuf.data(), MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, icomm, &stat);
            assert(0==ierr);

            int count;
            MPI_Get_count(&stat, MPI_CHAR, &count); // size of received data
            string const chbuf(recvbuf,0,count);    // C++
            cout << "   " << stat.MPI_SOURCE << " runs on " << chbuf;
            cout << " (length: " << count << " )" << endl;
            //  stat.Get_error()	// Error code
        }
    }
    else {
        int dest = 0;
        ierr = MPI_Send(name.data(), name.size(), MPI_CHAR, dest, myrank, icomm);
        assert(0==ierr);
    }
    return;
}
