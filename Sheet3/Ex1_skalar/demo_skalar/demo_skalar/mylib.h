#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <cassert>
#include <iomanip>             // setw()
#include <iostream>
#include <vector>

/** 	Inner product
	@param[in] x	vector
	@param[in] y	vector
	@return 	    resulting Euclidian inner product <x,y>
*/
double scalar(std::vector<double> const &x, std::vector<double> const &y);


//#########################################################################
// Defined by us
/** 	Inner product parallized without for loop
	@param[in] x	vector
	@param[in] y	vector
	@return 	    resulting Euclidian inner product <x,y>
*/
double scalarNofor(std::vector<double> const &x, std::vector<double> const &y);

/** 	Create vectors with iota and append them orderly
	@param[in] n	the length of each vector
	@return 	    resulting a large appended vector
*/
std::vector<int> reduction_vec_append_manually(int n);

/** 	Create vectors with iota and append them non-orderly
	@param[in] n	the length of each vector
	@return 	    resulting a large ordered appended vector
*/
std::vector<int> reduction_vec_append(int n);

/** 	 Append vector omp_in to omp_out .
	@param[in] omp_in	vector
	@param[in] omp_out	vector
	@return 	     vector with all omp_in appended 
*/
#pragma omp declare reduction(OurAppend : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) initializer (omp_priv(omp_orig))
//#########################################################################



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

/*
template <class T>
#pragma omp declare reduction(VecAdd : std::vector<T>  : omp_out += omp_in) initializer (omp_priv(omp_orig))
*/
// ------------------------------------------------------------


/** 	Test for vector reduction.
 * 
 * The thread-private vectors of size @p n are initialized via @f$v_k^{tID}=tID+k@f$.
 * Afterwards these vectors are accumulated, i.e., 
 * @f$v_k= \sum_{tID=0}^{numThreads} v_k^{tID}@f$.
 * 
 * 	@param[in] n  size of global/private vector
 * 	@return  resulting global vector.
*/
std::vector<int> reduction_vec(int n);



/** 	Output of a vector.
	@param[in,out] s	output stream
	@param[in]     x	vector
	@return 	   modified output stream
*/
template <class T>
std::ostream &operator<<(std::ostream &s, std::vector<T> const &x)
{
    for (auto const &v : x)  s << std::setw(4) << v << "  ";
    return s;
}

#endif
