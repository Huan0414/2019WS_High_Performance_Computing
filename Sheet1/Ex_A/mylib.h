#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <vector>

/** 	Inner product
	@param[in] x	vector
	@param[in] y	vector
	@return 	    resulting Euclidian inner product <x,y>
*/
double scalar(std::vector<double> const &x, std::vector<double> const &y);

/** 	L_2 Norm of a vector
	@param[in] x	vector
	@return 	    resulting Euclidian norm <x,y>
*/
double norm(std::vector<double> const &x);

vector<double> calMean(vector<double> const &x);


#endif
