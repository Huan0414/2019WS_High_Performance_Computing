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

    //const auto procx = static_cast<int>(sqrt(numprocs + 0.0));
    //const int  procy = procx;

    //if (procy * procx != numprocs)
    //{
        //cout << "\n Wrong number of processors !\n \n";
    //}
    //else
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
        //int nrefine = 1;
        //if (argc>1)  nrefine = atoi(argv[1]);

        ////Mesh const mesh_c("square_4.txt");
        //ParMesh const mesh_c("square");
        ////mesh_c.Debug();
        ////mesh_c.DebugEdgeBased();

        ////RefinedMesh mesh(mesh_c);                  // OK, works
        //////mesh.Debug();
        ////mesh.RefineAllElements(nrefine);           // OK, works

        //gMesh_Hierarchy ggm(mesh_c,nrefine);
        //const Mesh& mesh=ggm.finest();
        ////const Mesh& mesh=mesh_c;
        ParMesh const mesh("square_bb");

        //mesh.Debug();
        //mesh.DebugEdgeBased();
        //assert(mesh.CheckInterfaceAveraging());


        FEM_Matrix SK(mesh);                   // CRS matrix
        //SK.Debug();

        vector<double> uv(SK.Nrows(),0.0);     // temperature
        vector<double> fv(SK.Nrows(),0.0);     // r.h.s.

        mesh.VecAccu(uv);                      // Check MPI comm.

        SK.CalculateLaplace(fv);
        //if (0==myrank) mesh.Visualize(fv);
        //cout << mesh << endl;
        //SK.Debug();
        //
        mesh.SetValues(uv, [](double x, double y) -> double {
            //return x * x * std::sin(2.5 * M_PI * y);
            //return x * y;
            return 1.0;
            } );
        
        SK.ApplyDirichletBC(uv,fv);

        double tstart = MPI_Wtime();                  // Wall clock

        JacobiSolve(SK, fv, uv );          // solve the system of equations
        //JacobiSolve(mesh, SK, fv, uv );          // MPI: solve the system of equations

        double t1 =  MPI_Wtime() - tstart;             // Wall clock
        cout << "JacobiSolve: timing in sec. : " << t1 << endl;

        mesh.DebugValueAtCoords(uv, 3.0, 3.045,1e-2);


        if (!mesh.IsVectorConsistent(uv))
        {
            auto idx = mesh.IndicesOfInconistentData(uv);
            //assert(0==idx.size());
        }

        //mesh.Write_ascii_matlab("uv.txt", uv);
        //if (2==myrank || (1==numprocs && 0==myrank) ) mesh.Mesh::Visualize(uv);  // Visualize only one subdomain
        mesh.Visualize(uv);       // Visualize all subdomains
    }

    MPI_Finalize();
    return 0;
}
}


        //-----
        //vector<double> tu(SK.Nrows(),1.0);
        //vector<double> tf(SK.Nrows(),-9876.0);
        //SK.Mult(tf,tu);
        //if (2==myrank) mesh.Visualize(tf);
        //mesh.VecAccu(tf);
        //if (2==myrank) mesh.Visualize(tf);
        //vector<double> tf(SK.Nrows(),0.0);
        //mesh.SetDirchletValues(tf, [](double x, double y) -> double {
            ////return x * x * std::sin(2.5 * M_PI * y);
            //return x * y;
            //} );
        //if (3==myrank) mesh.Visualize(tf);
        //mesh.VecAccu(tf);
        //if (3==myrank) mesh.Visualize(tf);

        //vector<double> dd(SK.Nrows());           // matrix diagonal
        //SK.GetDiag(dd);                          //  dd := diag(K)
        //mesh.VecAccu(dd);                        //  MPI comm.  // take care for 1.0 for Dirichlet nodes
        //vector<double> tu(SK.Nrows(),1.0);
        //mesh.VecAccu(tu);
        //vddiv(tu,tu,dd);
        ////DebugVector(dd,"DD: ");{int ijk; cin >> ijk;}
        //if (3==myrank) mesh.Visualize(tu);
        //-----



