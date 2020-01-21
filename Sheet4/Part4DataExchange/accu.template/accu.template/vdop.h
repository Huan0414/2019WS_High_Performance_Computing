#ifndef VDOP_FILE
#define VDOP_FILE
#include <iostream>
#include <mpi.h>                  // MPI
#include <string>
#include <vector>

/** @brief  Element-wise vector divison x_k = y_k/z_k.
 *
 * @param[out] x  target vector
 * @param[in]  y  source vector
 * @param[in]  z  source vector
 *
 */
void vddiv(std::vector<double> &x, std::vector<double> const &y,
           std::vector<double> const &z);

/** @brief  Element-wise daxpy operation x(k) = y(k) + alpha*z(k).
 *
 * @param[out] x  target vector
 * @param[in]  y  source vector
 * @param[in]  alpha  scalar
 * @param[in]  z  source vector
 *
 */
void vdaxpy(std::vector<double> &x, std::vector<double> const &y,
            double alpha, std::vector<double> const &z );


/** @brief  Calculates the Euclidean inner product of two vectors.
 *
 * @param[in]  x vector
 * @param[in]  y vector
 * @return Euclidean inner product @f$\langle x,y \rangle@f$
 *
 */
double dscapr(std::vector<double> const &x, std::vector<double> const &y);


inline
double L2_scapr(std::vector<double> const &x, std::vector<double> const &y)
{
    return dscapr(x, y) / x.size();
}


/** 	Parallel inner product
  @param[in] x	    vector
  @param[in] y	    vector
  @param[in] icomm  MPI communicator
  @return 	        resulting Euclidian inner product <x,y>
*/
double par_scalar(std::vector<double> const &x, std::vector<double> const &y,
                   MPI_Comm const& icomm=MPI_COMM_WORLD);



/*  ReadId : Input and broadcast of an integer */
inline
int ReadIn(std::string const &ss = std::string(), MPI_Comm const &icomm = MPI_COMM_WORLD)
{
    MPI_Barrier(icomm);
    int myrank;		                    /*  my rank number */
    MPI_Comm_rank(icomm, &myrank);
    int id;

    if (myrank == 0) {
        std::cout << "\n\n  " << ss << " :  Which process do you want to debug ? \n";
        std::cin >> id;
    }
    MPI_Bcast(&id, 1, MPI_INT, 0, icomm);

    return id;
}

/**
 * Print entries of a vector.
 * @param[in] v	    vector values
*/
//void DebugVector(std::vector<double> const &v);
template <class T>
void DebugVector(std::vector<T> const &v, std::string const &ss = std::string(), MPI_Comm const &icomm = MPI_COMM_WORLD)
{
    MPI_Barrier(icomm);
    int numprocs;		                /*  # processes    */
    MPI_Comm_size(icomm, &numprocs);
    int myrank;		                    /*  my rank number */
    MPI_Comm_rank(icomm, &myrank);

    int readid = ReadIn(ss);				/* Read  readid    */

    while ( (0 <= readid) && (readid < numprocs) ) {
        if (myrank == readid) {
            std::cout << "\n\n process " << readid;
            std::cout << "\n ....  " << ss << " (nnode = " << v.size() << ")\n";
            for (size_t j = 0; j < v.size(); ++j) {
                std::cout.setf(std::ios::right, std::ios::adjustfield);
                std::cout << v[j] << "   ";
            }
            std::cout << std::endl;
            fflush(stdout);
        }

        readid = ReadIn(ss, icomm);			/* Read  readid    */
    }
    MPI_Barrier(icomm);
    return;
}

/** @brief  Compares an STL vector with POD vector.
 *
 * The accuracy criteria @f$ |x_k-y_k| < \varepsilon \left({1+0.5(|x_k|+|y_k|)}\right) @f$
 * follows the book by
 * <a href="https://www.springer.com/la/book/9783319446592">Stoyan/Baran</a>, p.8.
 *
 * @param[in]  x    STL vector
 * @param[in]  n    length of POD vector
 * @param[in]  y    POD vector
 * @param[in]  eps  relative accuracy criteria (default := 0.0).
 * @return true iff pairwise vector elements are relatively close to each other.
 *
 */
bool CompareVectors(std::vector<double> const &x, int n, double const y[], double const eps = 0.0);


/** Output operator for vector
 * 	@param[in,out] s	output stream, e.g. @p cout
 *  @param[in]     v    vector
 *
 *	@return    output stream
*/
template <class T>
std::ostream &operator<<(std::ostream &s, std::vector<T> const &v)
{
    for (auto vp : v) {
        s << vp << "  ";
    }
    return s;
}

/** Exchanges equal size partions of vector @p xin with all MPI processes.
 * The received data are return in vector @p yout .
 * 
 *  @param[in]  xin   input vector
 *  @param[out] yout  output vector
 *  @param[in]  icomm MPI communicator
 *
*/
void ExchangeAll(std::vector<double> const &xin, std::vector<double> &yout, MPI_Comm const &icomm = MPI_COMM_WORLD);

/** Exchanges equal size partions of vector @p xin with all MPI processes.
 * The received data are return in vector @p xin  .
 * 
 *  @param[in,out]  xin   input/output vector
 *  @param[in]  icomm MPI communicator
 *
*/
void ExchangeAllInPlace(std::vector<double> &xin, MPI_Comm const &icomm = MPI_COMM_WORLD);



#endif
