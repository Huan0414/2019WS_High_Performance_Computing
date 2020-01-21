#include "greetings.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <mpi.h>                                  // MPI
#include <string>
using namespace std;

// see  http://www.open-mpi.org/doc/current
// for details on MPI functions

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

            ierr = MPI_Recv(chbuf, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, icomm, &stat);
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
