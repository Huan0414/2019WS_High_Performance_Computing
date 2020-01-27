//	general header for all functions in directory

#ifndef GREETINGS_FILE
#define GREETINGS_FILE

#include <mpi.h>
#include <vector>

/**	Each process finds out its host, sends this information
	to root process 0 which prints this information for each process.
	@param[in]	icomm	the MPI process group that is used.
*/

void greetings(MPI_Comm const &icomm);
void greetings_cpp(MPI_Comm const &icomm);
void DebugVector(std::vector<double> const &xin, MPI_Comm const &icomm );
void PrintVectorsInOrder(std::vector<double> const &xin, MPI_Comm const &icomm );
void scalarProduct(double &sum, std::vector<double> const &v1, std::vector<double> const &v2, MPI_Comm const &icomm);
void par_scalar(double &innerProduct, std::vector<double> const &vec1, std::vector<double> const &vec2, MPI_Comm const &icomm);
void MinMaxWithoutReduc(std::vector<double> &localVec, MPI_Comm const icomm);

#endif
