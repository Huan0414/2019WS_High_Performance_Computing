//	general header for all functions in directory

#ifndef GREETINGS_FILE
#define GREETINGS_FILE

#include <mpi.h>

/**	Each process finds out its host, sends this information
	to root process 0 which prints this information for each process.
	@param[in]	icomm	the MPI process group that is used.
*/

void greetings(MPI_Comm const &icomm);
void greetings_cpp(MPI_Comm const &icomm);

#endif
