#include "log_fft.h"
#include <stdlib.h>
//  utils is included from log_fft

//  #define _DEBUG_PRINT_
//  #define _SUPER_DEBUG_

void test_log_fft(double* input_vec, double* output_vec, int n)
{
    complex_sl_t* buf = (complex_sl_t*) malloc(n * sizeof(complex_sl_t));

    for(int i1=0; i1<n; i1++){
        buf[i1].real = quantizer(input_vec[i1]);
        buf[i1].imag = quantizer(0);
    }

#ifdef _DEBUG_PRINT_
    print_array_cplx_sl("Data: ", buf, n);
#endif

    sl_fft(buf, n);

#ifdef _DEBUG_PRINT_
    print_array_cplx_sl("FFT : ", buf, n);
#endif

    for(int i=0; i<n; i++)  output_vec[i] = inverse(buf[i].real);
    
    free(buf);
}

void simple_test(){
    double input_vec[8] = {0.250, -0.178, 0.511, 0.006, 0.557, -0.883, -0.374, -0.295};
    double output_vec1[8] = {1.0,1.0,1.0,1.0,0,0,0,0};
    double output_vec2[8] = {1.0,1.0,1.0,1.0,0,0,0,0};
    test_log_fft(input_vec, output_vec2, 8);
}

int main(){
#ifdef HW_COSIM
    simple_test();
#endif
	/*
	 *
	for(int i=1; i<12; i++){
        test_n_point_fft(1 << i);
    }
    */
    
    return 0;
}
