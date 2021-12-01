#include "log_fft.c"
#include <stdlib.h>
//  utils is included from log_fft

//  #define _DEBUG_PRINT_
//  #define _SUPER_DEBUG_

void test_log_fft(double* input_vec, double* output_vec)
{
    complex_sl_t buf[FFT_POINT];

    for(int i1=0; i1<FFT_POINT; i1++){
        
        sig_log_t real;
        sig_log_t imag;
        quantizer(input_vec[i1], real);
        quantizer(0, imag);

        copy_two_comp_to_cplx(buf[i1], real, imag);
    }

#ifdef _DEBUG_PRINT_
    print_array_cplx_sl("Data: ", buf, n);
#endif

    sl_fft(buf);

#ifdef _DEBUG_PRINT_
    print_array_cplx_sl("FFT : ", buf, n);
#endif

    for(int i=0; i<FFT_POINT; i++){
        sig_log_t real;
        sig_log_t imag;
        copy_cplx_to_two_comp(real, imag, buf[i]);
        output_vec[i] = inverse(imag);
    }
    //  free(buf);
}

void simple_test(){
    double input_vec[FFT_POINT];
    for(int i=0; i<FFT_POINT; i++)  input_vec[i] = randfrom(0,1);
    double output_vec1[FFT_POINT];
    test_log_fft(input_vec, output_vec1);
    for(int i=0; i<FFT_POINT; i++)  printf("%g, ", output_vec1[i]);
    printf("What's up\n");
}

int main(){
//#ifdef HW_COSIM
    simple_test();
//#endif
	/*
	 *
	for(int i=1; i<12; i++){
        test_n_point_fft(1 << i);
    }
    */
    printf("What's up");
    return 0;
}
