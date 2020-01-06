#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <cassert>
#include <iomanip>             // setw()
#include <iostream>
#include <vector>

/**	Calculate minmum, maximum, arithmetic mean, geometric mean, harmonic mean and standard deviation
	@param[in] n 
	@param[in] pair_cont vector, initialized as 0
	@retun 		resulting new pair_cont
*/
void count_goldbach(const int n, std::vector<int> &pair_cont, std::vector<int> const &primes);
void count_goldbachTwo(const int n, std::vector<int> &pair_cont, std::vector<int> const &primes);
void count_goldbachCheck(const int n, std::vector<int> &pair_cont, std::vector<int> const &primes);

std::vector<int> get_primes(int max);


/** 	 Vector @p b adds its elements to vector @p a .
	@param[in] a	vector
	@param[in] b	vector
	@return 	     a+=b componentwise
*/
template<class T>
std::vector<T> &operator+=(std::vector<T> &a, std::vector<T> const &b)
{
    assert(a.size()==b.size());
    for (size_t k = 0; k < a.size(); ++k) {
        a[k] += b[k];
    }
    return a;
}

// Declare the reduction operation in OpenMP for an STL-vector
//   omp_out += omp_in  requires operator+=(vector<int> &, vector<int> const &) from above
//   !!  CLANG_  does not compile the pragma declare
// ------------------------------------------------------------
// https://scc.ustc.edu.cn/zlsc/tc4600/intel/2016.0.109/compiler_c/common/core/GUID-7312910C-D175-4544-99C5-29C12D980744.htm
// https://gist.github.com/eruffaldi/7180bdec4c8c9a11f019dd0ba9a2d68c
// https://stackoverflow.com/questions/29633531/user-defined-reduction-on-vector-of-varying-size
//  see also p.74ff in  https://www.fz-juelich.de/ias/jsc/EN/AboutUs/Staff/Hagemeier_A/docs-parallel-programming/OpenMP-Slides.pdf
#pragma omp declare reduction(VecAdd : std::vector<int>  : omp_out += omp_in) initializer (omp_priv(omp_orig))
//#pragma omp declare reduction(VecAdd : std::vector<int>  : omp_out += omp_in) initializer (omp_priv=omp_orig) // compiles with CLANG_ but fails at runtime



#endif
