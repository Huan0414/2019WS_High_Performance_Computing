#ifndef FILE_MYLIB               // header guarding
#define FILE_MYLIB

// Function with  o n e  return value
/** Calculates Kelvin from degree Celsius.
*
*  @param[in]  grad_C Temperature in degree Celsius
*
*  @return Temperature in Kelvin
*  @warning No check for non-physical values (<0 Kelvin)
*/
float c2k(float const grad_C);            // declaration of function


// Function with  m u l t i p le   return values
/** Calculates Kelvin and degree Fahrenheit from degree Celsius.
*
*  @param[in]  grad_C Temperature in degree Celsius
*  @param[out] grad_K Temperature in Kelvin
*  @param[out] grad_F Temperature in degree Fahrenheit
*
*  @warning Code stops for non-physical values (<0 Kelvin).
*           Compiler option -NDEBUG  switches off this check.
*/
void c2kf(float const grad_C, float& grad_K, float &grad_F); // declaration of function

#endif
