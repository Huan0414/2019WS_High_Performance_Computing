#ifndef HEADER_CHECK_ENV
#define HEADER_CHECK_ENV

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <unordered_map>

//#####################################
// G.Haase
// See https://sourceforge.net/p/predef/wiki/Compilers/
//     http://www.cplusplus.com/doc/tutorial/preprocessor/
//#####################################
/** 	Checks for compilers, its versions, threads etc.
 * 
	@param[in] argc	number of command line arguemnts
	@param[in] argv	command line arguemnts as array of C-strings
*/
template <class T>
void check_env(T argc, char const *argv[])
{
    std::cout << "\n#######################################################################\n";
    std::cout << "Code    :";
    for (T k = 0; k < argc; ++k) std::cout << "  " << argv[k];
    std::cout << std::endl;

    // compiler:      https://sourceforge.net/p/predef/wiki/Compilers/
    std::cout <<    "Compiler:  ";
#if defined __INTEL_COMPILER
#pragma message(" ##########  INTEL  ###############")
    std::cout << "INTEL " << __INTEL_COMPILER;
    // Ignore warnings for #pragma acc   unrecognice
#pragma warning disable 161
    // Ignore warnings for #pragma omp   unrecognice
#pragma warning disable 3180

#elif defined __PGI
#pragma message(" ##########  PGI    ###############")
    std::cout << "PGI " << __PGIC__ << "." << __PGIC_MINOR__ << "." << __PGIC_PATCHLEVEL__;
#elif defined  __clang__
#pragma message(" ##########  CLANG    ###############")
    std::cout << "CLANG " << __clang_major__ << "." << __clang_minor__ << "."; // << __clang_patchlevel__;
#elif defined __GNUC__
#pragma message(" ##########  Gnu    ###############")
    std::cout << "Gnu " <<  __GNUC__  << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
#pragma message(" ##########  unknown Compiler   ###############")
    std::cout << "unknown";
#endif
    std::cout << "  C++ standard: " << __cplusplus << std::endl;

    // Parallel environments
    std::cout <<    "Parallel:  ";
#if defined MPI_VERSION
#pragma message(" ##########  MPI    ###############")
#ifdef OPEN_MPI
    std::cout << "OpenMPI ";
#else
    std::cout << "MPI ";
#endif
    std::cout << MPI_VERSION << "." << MPI_SUBVERSION << "   ";
#endif

#ifdef _OPENMP
//https://www.openmp.org/specifications/
//https://stackoverflow.com/questions/1304363/how-to-check-the-version-of-openmp-on-linux
    std::unordered_map<unsigned, std::string> map{
        {200505, "2.5"}, {200805, "3.0"}, {201107, "3.1"}, {201307, "4.0"}, {201511, "4.5"}, {201811, "5.0"}};
#pragma message(" ##########  OPENMP    ###############")
    //std::cout << _OPENMP;
    std::cout << "OpenMP ";
    try {
        std::cout << map.at(_OPENMP);
    }
    catch (...) {
        std::cout << _OPENMP;
    }
    #pragma omp parallel
    {
        #pragma omp master
        {
            const int nn = omp_get_num_threads();          // OpenMP
            std::cout << " ---> " <<  nn << " Threads   ";
        }
        #pragma omp barrier
    }

#endif
#ifdef _OPENACC
#pragma message(" ##########  OPENACC    ###############")
    std::cout << "OpenACC   ";
#endif
    std::cout << std::endl;
    std::cout << "Date    :  " << __DATE__ << "  " << __TIME__;
    std::cout << "\n#######################################################################\n";
}
// HG

#endif
