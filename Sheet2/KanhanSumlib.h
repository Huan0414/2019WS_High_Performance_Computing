#ifndef FILE_MYLIB
#define FILE_MYLIB
#include <vector>
#include <cmath>
#include <limits>
#define M_PI    3.14159265358979323846 /*pi*/

double scalar(std::vector<double> const &a, long long int const N);
double Kahan_scalar(std::vector<double> const &a, long long int const N);

#endif
