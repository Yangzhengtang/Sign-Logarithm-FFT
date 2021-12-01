#include "log_fft.c"
#include <stdlib.h>
//  utils is included from log_fft

//  #define _DEBUG_PRINT_
//  #define _SUPER_DEBUG_

void test_log_fft(double* input_vec, double* output_vec, int n)
{
    complex_sl_t buf[8];

    for(int i1=0; i1<n; i1++){
        //  buf[i1].real = quantizer(input_vec[i1]);
        //  buf[i1].imag = quantizer(0);
        
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

    for(int i=0; i<n; i++){
        sig_log_t real;
        sig_log_t imag;
        copy_cplx_to_two_comp(real, imag, buf[i]);
        output_vec[i] = inverse(imag);
    }
    //  free(buf);
}

void simple_test(){
    double input_vec[8] = {1, 1, 1, 1, 0, 0, 0, 0};
    double output_vec1[8] = {1.0,1.0,1.0,1.0,0,0,0,0};
    test_log_fft(input_vec, output_vec1, 8);
    for(int i=0; i<9; i++)  printf("%g, ", output_vec1[i]);
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
