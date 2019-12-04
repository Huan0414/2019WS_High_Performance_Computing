#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <vector>

/** 	Matrix-vector product
 * 	@param[in] A	dense matrix (2D access)
 *  @param[in] u	vector
 *
 *	@return    resulting vector
 *  @warning   Non continuous memory access for matrix.
 *  @warning   It is assumed that rows in @p A have the same number of elements (no check!).
*/
std::vector<double> MatVec(std::vector<std::vector<double>> const &A, std::vector<double> const &u);


/** 	Matrix-vector product
 * 	@param[in] A	dense matrix (1D access)
 *  @param[in] u	vector
 *
 *	@return    resulting vector
*/
std::vector<double> MatVec(std::vector<double> const &A, std::vector<double> const &u);



/** 	Operator overloading: Matrix-vector product with the multiplication operator
 * 	@param[in] A	dense matrix (2D access)
 *  @param[in] u	vector
 *
 *	@return    resulting vector
 *  @warning   Non continuous memory access for matrix.
 *  @warning   It is assumed that rows in @p A have the same number of elements (no check!).
*/
inline
std::vector<double> operator*(std::vector<std::vector<double>> const &A, std::vector<double> const &u)
{
    return MatVec(A,u);
}

/** 	Operator overloading: Matrix-vector product with the multiplication operator
 * 	@param[in] A	dense matrix (1D access)
 *  @param[in] u	vector
 *
 *	@return    resulting vector
*/
inline
std::vector<double> operator*(std::vector<double> const &A, std::vector<double> const &u)
{
    return MatVec(A,u);
}


#endif
