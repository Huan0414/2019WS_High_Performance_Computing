//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "geom.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"
#include "vdop.h"

#include <cmath>
#include <iostream>
#include <omp.h>
#include <ctime>
#include <sstream>
#include <string>

using namespace std;


int main(int , char ** )
{
    const int numprocs = 1;
    const int myrank   = 0;

    if (myrank == 0)
    {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }

    const auto procx = static_cast<int>(sqrt(numprocs + 0.0));
    const int  procy = procx;

    if (procy * procx != numprocs)
    {
        cout << "\n Wrong number of processors !\n \n";
    }
    else
    {
// #####################################################################
//      Here starts the real code
// #####################################################################
        //bool ScaleUp = !true;
        int nx, ny, NXglob, NYglob; /* number of local intervals on (xl,xr)=:nx, (yb,yt)=:ny */
        //nx = 1024;
        //ny = 1024;
        nx = 1000;
        ny = 1000;
        NXglob = nx * procx;
        NYglob = ny * procy;
        cout << "Intervalls: " << NXglob << " x " << NYglob << endl;

// ##################### STL ###########################################
{
        ////########################################################################
    int nthreads;                                  // OpenMP
    #pragma omp parallel default(none) shared(cout,nthreads)
    {
        //int const th_id  = omp_get_thread_num();   // OpenMP
        int const nthrds = omp_get_num_threads();  // OpenMP
        //stringstream ss;
        //ss << "C++: Hello World from thread " << th_id << " / " << nthrds << endl;
        //#pragma omp critical
        //{
            //cout << ss.str();                      // output to a shared ressource
        //}
        //#pragma omp master
        nthreads = nthrds;                         // transfer nn to to master thread
    }
    cout << "   " << nthreads << "   threads have been started." << endl;

        
        Mesh_2d_3_square const mesh(nx, ny);
        //mesh.Debug();

        CRS_Matrix SK(mesh);                   // CRS matrix
        //SK.Debug();
        vector<double> uv(SK.Nrows(),0.0);     // temperature
        vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        SK.CalculateLaplace(fv);
        //SK.Debug();

        //mesh.SetU(uv);         // deprecated
        //mesh.SetF(fv);         // deprecated
        // Two ways to initialize the vector
        //mesh.SetValues(uv,f_zero);             // functional
        mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function

        SK.ApplyDirichletBC(uv,fv);
        //SK.Compare2Old(nnode, id, ik, sk);
        //SK.Debug();

        //double tstart = clock();                        // timing
        double tstart = omp_get_wtime();

        JacobiSolve(SK, fv, uv );          // solve the system of equations

        //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
        double t1 = omp_get_wtime() - tstart;
        cout << "JacobiSolve: timing in sec. : " << t1 << endl;

        //CompareVectors(uv, nnode, u, 1e-6);    // Check correctness

        //mesh.SaveVectorP("t.dat", uv);
        //mesh.Visualize(uv);

}
// ##################### STL ###########################################
{
        //Mesh_2d_3_matlab const mesh("square_tiny.txt");
        Mesh_2d_3_matlab const mesh("square_100.txt");
        //Mesh_2d_3_matlab const mesh("L_shape.txt");
        //mesh.Debug();

        CRS_Matrix SK(mesh);                   // CRS matrix
        //SK.Debug();

        vector<double> uv(SK.Nrows(),0.0);     // temperature
        vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        SK.CalculateLaplace(fv);
        //SK.Debug();

        //mesh.SetU(uv);         // deprecated
        // Two ways to initialize the vector
        //mesh.SetValues(uv,f_zero);             // user function
        mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function

        SK.ApplyDirichletBC(uv,fv);
        //SK.Compare2Old(nnode, id, ik, sk);
        //SK.Debug();

        //double tstart = clock();                        // timing
		double tstart = omp_get_wtime();
        JacobiSolve(SK, fv, uv );          // solve the system of equations

        //double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
        double t1 = omp_get_wtime() - tstart;
        cout << "JacobiSolve: timing in sec. : " << t1 << endl;

        //mesh.Write_ascii_matlab("uv.txt", uv);
        //mesh.Visualize(uv);
    }
    return 0;
}
}


