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

       assert(4 == np);                  // example is only provided for 4 MPI processes
   }
// #####################################################################
// ---- Read the f.e. mesh and the mapping of elements to MPI processes
    //Mesh const mesh_c("square_4.txt");    //    Files square_4.m and square_4_sd.txt  are needed
    ParMesh const mesh("square",icomm);
    
    int const numprocs = mesh.NumProcs();
    int const myrank   = mesh.MyRank();
    if ( 0 == myrank ) {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }
        
    int const check_rank=0;               // choose the MPI process you would like to check the mesh
    //if ( check_rank == myrank ) mesh.Debug();
    //if ( check_rank == myrank ) mesh.DebugEdgeBased();
    
// ---- allocate local vectors and check skalar product and vector accumulation
    vector<double> xl(mesh.Nnodes(), -1.0);
    mesh.SetValues(xl, [](double x, double y) -> double {return x * x * std::sin(2.5 * M_PI * y);} );
    
    //if (check_rank==myrank) mesh.Visualize(xl);
    
    for (size_t k=0; k<xl.size(); ++k)
    {
        xl[k] = 1.0;
    }
    double ss = mesh.dscapr(xl,xl);
    cout << myrank << " : scalar : " << ss << endl;
    
    mesh.VecAccu(xl);
    
    if (check_rank==myrank) mesh.Visualize(xl);

    MPI_Finalize();
    return 0;
}


