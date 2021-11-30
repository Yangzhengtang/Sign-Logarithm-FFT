//  Simple test
#include "log_fft.h"

int simple_test(){
    fixed_point_t a;
    fixed_point_t b;
    fixed_point_t ret;
    a[0] = 1;
    a[1] = 1;
    b[1] = 2;
    b[1] = 2;
    multi(a, b, ret);
}

int main(){
#ifdef HW_COSIM
    simple_test();
#endif
    return 0;
}