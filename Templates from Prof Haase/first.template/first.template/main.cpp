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
    numprocs = 1;    // delete this line when uncommenting the next line
    MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);

	if (0==myrank) {
        cout << "\n There are " << numprocs << " processes running.\n \n";
    }

    greetings(icomm);
    greetings_cpp(icomm);

    if (0==myrank) cout << endl;

    MPI_Finalize();

    return 0;
}


