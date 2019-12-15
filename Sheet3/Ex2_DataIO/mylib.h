#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <vector>
#include <cassert>
#include <algorithm>

/**	Calculate minmum, maximum, arithmetic mean, geometric mean, harmonic mean and standard deviation
	@param[in] a vector
	@param[in] b vector
	@retun 		resulting b(6)
*/
//void cal_min_max_mean(std::vector<double> const &a, std::vector<double> &b);
std::vector<double> cal_min_max_mean(std::vector<double> const &a);

template<class T>
std::vector<T> &operator+=(std::vector<T> &a, std::vector<T> const &b) 
	{
	assert((a.size() == 6)&&(b.size()==6));
	a[0] = std::min(a[0],b[0]);
	a[1] = std::max(a[1],b[1]);
	a[2] += b[2];
	a[3] *= b[3];
	a[4] += b[4];
	a[5] += b[5];
	return a;
	}

#pragma omp declare reduction (ReductionForMeansCalculation : std::vector<double> : omp_out += omp_in) initializer (omp_priv(omp_orig))

#endif
