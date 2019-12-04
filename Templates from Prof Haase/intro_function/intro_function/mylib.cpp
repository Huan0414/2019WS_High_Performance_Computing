#include "mylib.h"
#include <cassert>          // assert

// Function with  o n e  return value
float c2k(float const grad_C)             // definition of function
{
    float grad_K;
    grad_K = grad_C + 273.15f;            // D o n ' t   d o !   273,15
    return grad_K;
}


// Function with  m u l t i p le   return values
void c2kf(float const grad_C, float& grad_K, float &grad_F)  // definition of function
{
    assert( grad_C >= -273.15f );         // Stop in case of non-physical values
    grad_K = c2k(grad_C);
    grad_F = 9/5.0f*grad_C + 32;          // D o n ' t   d o !   9/5
    return;
}
