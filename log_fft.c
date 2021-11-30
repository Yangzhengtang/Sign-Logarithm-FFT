#ifndef _LOGFFT_INCL_
#define _LOGFFT_INCL_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "log_fft.h"
#include <time.h>

void multi(fixed_point_t a, fixed_point_t b, fixed_point_t ret){
    ret[0] = a[0] * b[0];
    ret[1] = a[1] + b[1];
    return;
}

void multi_cplx(complex_sl_t a, complex_sl_t b, complex_sl_t ret){
    multi(a[0], b[0], ret[0]);
}

#endif
