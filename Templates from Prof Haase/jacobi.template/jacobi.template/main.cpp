//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "geom.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "par_geom.h"
#include "userset.h"
#include "vdop.h"

#include <cmath>
#include <iostream>
#include <mpi.h>            // MPI
#include <omp.h>            // OpenMP
using namespace std;


int main(int argc , char **argv )
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
    ParMesh const mesh("square");

    int const numprocs = mesh.NumProcs();
    int const myrank   = mesh.MyRank();
    if ( 0 == myrank )
    {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }

    int const check_rank = 0;             // choose the MPI process you would like to check the mesh
    //if ( check_rank == myrank ) mesh.Debug();
    //if ( check_rank == myrank ) mesh.DebugEdgeBased();


    FEM_Matrix SK(mesh);                   // CRS matrix
    //SK.Debug();

    vector<double> uv(SK.Nrows(), 0.0);    // temperature
    vector<double> fv(SK.Nrows(), 0.0);    // r.h.s.

    mesh.VecAccu(uv);                      // Check MPI comm.

    SK.CalculateLaplace(fv);
    //SK.Debug();
    //
    mesh.SetValues(uv, [](double x, double y) -> double
    {
        return x *x * std::sin(2.5 * M_PI * y);
    } );

    ////mesh.SetU(uv);         // deprecated
    //// Two ways to initialize the vector
    ////mesh.SetValues(uv,f_zero);             // user function
    ////mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function
    ////mesh.SetValues(uv, [](double x, double y) -> double {return 5e-3*(x+1)*(y+1);} );  // lambda function
    ////
    //mesh.SetValues(uv, [](double x, double y) -> double {
    //return x * x * std::sin(2.5 * M_PI * y);
    //} );

    SK.ApplyDirichletBC(uv, fv);
    //SK.Debug();

    double tstart = MPI_Wtime();                  // Wall clock

    JacobiSolve(SK, fv, uv );          // solve the system of equations
    //JacobiSolve(mesh, SK, fv, uv );          // MPI: solve the system of equations

    double t1 =  MPI_Wtime() - tstart;             // Wall clock
    cout << "JacobiSolve: timing in sec. : " << t1 << endl;

    //if (check_rank == myrank) mesh.Visualize(uv);

    MPI_Finalize();
    return 0;
}


