//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "geom.h"
#include "par_geom.h"
#include "vdop.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>            // MPI
#include <omp.h>            // OpenMP
using namespace std;


int main(int argc, char **argv )
{
    MPI_Init(&argc, &argv);
    MPI_Comm const icomm(MPI_COMM_WORLD);
    omp_set_num_threads(1);                   // don't use OMP parallelization for a start
//
   {
       int np;
       MPI_Comm_size(icomm, &np);

       assert(6 == np);                  // example is only provided for 4 MPI processes
   }
// #####################################################################
// ---- Read the f.e. mesh and the mapping of elements to MPI processes
    //Mesh const mesh_c("square_4.txt");    //    Files square_4.m and square_4_sd.txt  are needed
    //ParMesh const mesh("square",icomm);
    //Mesh const mesh_c("Rectangle_6.txt");    //    Files Rectangle_6.m and Rectangle_6_sd.txt  are needed
    ParMesh const mesh("Rectangle",icomm);
    
    int const numprocs = mesh.NumProcs();
    int const myrank   = mesh.MyRank();
    if ( 0 == myrank ) {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }
        
    int const check_rank=0;               // choose the MPI process you would like to check the mesh
    //if ( check_rank == myrank ) mesh.Debug();
    //if ( check_rank == myrank ) mesh.DebugEdgeBased();
    
// ---- allocate local vectors and check skalar product and vector accumulation
    //vector<double> xl(mesh.Nnodes(), 1.0);
    vector<double> xl(mesh.Nnodes(), myrank + 1.0);
    vector<double> xl2(mesh.Nnodes(), myrank + 2.0); 
    vector<int> xlInt(mesh.Nnodes(), myrank + 1);		 
       
    //mesh.SetValues(xl, [](double x, double y) -> double {return x * x * std::sin(2.5 * M_PI * y);} );
    //mesh.SetValues( xl, [](double x, double y) -> double {return x + y;} );
    
    //if (check_rank==myrank) mesh.Visualize(xl);
    
    //for (size_t k=0; k<xl.size(); ++k)
    //{
    //    xl[k] = 1.0;
    //}
    
//// ####################### Sheet4Ex9 check VecAccu() and GetCoords() ###     
        
    if(myrank == check_rank) {
	cout << "Sheet4Ex9 check VecAccu() and GetCoords()" << endl;
    cout << "Nnodes: " << xl.size() << endl;
    }			
    cout << "MyRank: " << myrank << ", Nnodes: " << xl.size() << endl;
    
    //cout << "Coordinates:" << mesh.GetCoords() << endl;
    //cout << "Size of mesh: " << (int)mesh.GetCoords().size() / 2 << endl;		// we don't understand why need coords
       
    //if (check_rank==myrank) mesh.Visualize(xl);					// check before VecAccu()
    
    //mesh.VecAccu(xl);
 
    //if (1==myrank) mesh.Visualize(xl);					// check after VecAccu()
																	// 2 on the intefaces, 4 on the center
    
    //double ss = mesh.dscapr(xl,xl);
	//cout << myrank << " : scalar : " << ss << endl;				//  should be 1108
																	////  overall 660 interior nodes, 48 interfaces nodes, 1 center node

//// ####################### Sheet4Ex10 VecAccuInt() with int ##############    
    
    //MPI_Barrier(icomm);
    
    //mesh.VecAccuInt(xlInt);    
	
	//int ssInt = mesh.dscaprInt(xlInt,xlInt);
    //cout << myrank << " : scalar : " << ssInt << endl;

    //if (check_rank==myrank) mesh.Visualize(xl);
    
    
//// ####################### Sheet4Ex11 GlobalNnodes() ###################   
    
    //MPI_Barrier(icomm);			
    
    //int GlobalNumberNodes = mesh.GlobalNnodes();					// MPI_Allreduce used
    //if(myrank == check_rank) 
    //{cout << "GlobalNodesNumber: " << GlobalNumberNodes << endl;} // should be 709, (try with 3 domains)													
	
//// ####################### Sheet4Ex12 VecAverage() ################## 
    
    //MPI_Barrier(icomm);
    
    mesh.VecAverage(xl);   
    if (check_rank==myrank) mesh.Visualize(xl);						// all be 1.0
    
    //mesh.VecAverage(xl2);  
    //if (check_rank==myrank) mesh.Visualize(xl2);					// all be 2.0
    
    MPI_Finalize();
    return 0;
}


