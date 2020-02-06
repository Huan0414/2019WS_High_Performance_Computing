//		MPI code in C++.
//		See [Gropp/Lusk/Skjellum, "Using MPI", p.33/41 etc.]
//		and  /opt/mpich/include/mpi2c++/comm.h  for details

#include "geom.h"
#include "par_geom.h"
#include "getmatrix.h"
#include "jacsolve.h"
#include "userset.h"
#include "vdop.h"

#include <cmath>
#include <iostream>
#include <mpi.h>            // MPI
#include <omp.h>            // OpenMP
using namespace std;


int main(int argc , char **argv )
{
    MPI_Init(&argc,&argv);
    MPI_Comm icomm(MPI_COMM_WORLD);
    int numprocs;
    MPI_Comm_size(icomm, &numprocs);
    int myrank;
    MPI_Comm_rank(icomm, &myrank);

    omp_set_num_threads(1);                   // don't use OMP parallelization for a start

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
        ////bool ScaleUp = !true;
        //int nx, ny, NXglob, NYglob; /* number of local intervals on (xl,xr)=:nx, (yb,yt)=:ny */
        ////nx = 1024;
        ////ny = 1024;
        //nx = 100;
        //ny = 100;
        //NXglob = nx * procx;
        //NYglob = ny * procy;
        //cout << "Intervalls: " << NXglob << " x " << NYglob << endl;

// ##################### STL ###########################################
{
        int nrefine = 1;
        if (argc>1)  nrefine = atoi(argv[1]);

        //Mesh const mesh_c("square_4.txt");
        ParMesh const mesh_c("square");
        //mesh_c.Debug();
        //mesh_c.DebugEdgeBased();

        //RefinedMesh mesh(mesh_c);                  // OK, works
        ////mesh.Debug();
        //mesh.RefineAllElements(nrefine);           // OK, works

        gMesh_Hierarchy ggm(mesh_c,nrefine);
        const Mesh& mesh=ggm.finest();
        //const Mesh& mesh=mesh_c;

        //mesh.Debug();
        //mesh.DebugEdgeBased();


        FEM_Matrix SK(mesh);                   // CRS matrix
        //SK.Debug();

        vector<double> uv(SK.Nrows(),0.0);     // temperature
        vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        SK.CalculateLaplace(fv);
        //SK.Debug();

        //mesh.SetU(uv);         // deprecated
        // Two ways to initialize the vector
        //mesh.SetValues(uv,f_zero);             // user function
        //mesh.SetValues(uv, [](double x, double y) -> double {return 0.0*x*y;} );  // lambda function
        //mesh.SetValues(uv, [](double x, double y) -> double {return 5e-3*(x+1)*(y+1);} );  // lambda function
        //
        mesh.SetValues(uv, [](double x, double y) -> double {
            return x * x * std::sin(2.5 * M_PI * y);
            } );


        SK.ApplyDirichletBC(uv,fv);
        //SK.Compare2Old(nnode, id, ik, sk);
        //SK.Debug();    MPI::Finalize();


        double tstart = clock();                        // timing

        JacobiSolve(SK, fv, uv );          // solve the system of equations

        double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
        cout << "JacobiSolve: timing in sec. : " << t1 << endl;

        //mesh.Write_ascii_matlab("uv.txt", uv);
        //if (0==myrank) mesh.Visualize(uv);
    }

    //{
        //int nrefine = 4;
        //if (argc>1)  nrefine = atoi(argv[1]);

        //Multigrid ggm(Mesh("square_4.txt"),nrefine);

        //ggm.DefineOperators();

        //cout << "\n#############  SOLVE   ###############\n";
        //int my_level=nrefine-1;

        ////double tstart = clock();                        // timing
        //double tstart = omp_get_wtime();                  // OpenMP

        ////ggm.JacobiSolve(my_level);
        ////ggm.MG_Step(my_level, 1, true, 1);
        //ggm.MG_Solve(2, 1e-6, 1);

        ////double t1 = (clock() - tstart) / CLOCKS_PER_SEC;// timing
        //double t1 = omp_get_wtime() - tstart;             // OpenMP
        //cout << "MgSolve: timing in sec. : " << t1 << "   for " << ggm.Ndofs()<< " dofs"<< endl;

        //const auto &ml=ggm.GetMesh(my_level);
        //const auto &sl=ggm.GetSolution(my_level);
        ////ml.Visualize(sl);

        //ml.Export_scicomp("level_"+to_string(my_level));
    //}

    MPI_Finalize();
    return 0;
}
}


