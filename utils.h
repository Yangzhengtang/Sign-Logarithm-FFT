/*
    Fixed - Floating convertion utilities,
    Reference: https://embeddedartistry.com/blog/2018/07/12/simple-fixed-point-conversion-in-c/
*/

#include <stdio.h>
#include <math.h>

/// Fixed-point Format: 3.13 (16-bit)
typedef long long fixed_point_t;
#define FIXED_POINT_FRACTIONAL_BITS 25
/// Converts 13.3 format -> double
double fixed_to_double(fixed_point_t input);

/// Converts double to 13.3 format
fixed_point_t double_to_fixed(double input);

inline double fixed_to_double(fixed_point_t input)
{
    return ((double)input / (double)(1 << FIXED_POINT_FRACTIONAL_BITS));
}

inline fixed_point_t double_to_fixed(double input)
{
    return (fixed_point_t)(round(input * (1 << FIXED_POINT_FRACTIONAL_BITS)));
}