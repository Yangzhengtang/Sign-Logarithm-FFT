#include "log_fft.h"
#include "fft.h"
//  utils is included from log_fft

//  #define _DEBUG_PRINT_
//  #define _SUPER_DEBUG_

void test_quantize_and_inverse(){
    double x = 0.31415926;
    char sig_x = 0;
    fixed_point_t  log_x = 0;
    sig_log_t sl_x = quantizer(x);

    printf("After quantizer, x: %f, sig_x: %d, log_x: %hu\n", x, sl_x.sig_a, sl_x.log_a);

    double new_x = inverse(sl_x);

    printf("new_x: %f\n", new_x);
}

void test_op(){
    double a = 7.5;
    double b = 7.5;
    sig_log_t sl_a = quantizer(a);
    sig_log_t sl_b = quantizer(b);

    printf("After quantizer, a: %f, sig_a: %d, log_a: %u\n", a, sl_a.sig_a, sl_a.log_a);
    printf("After quantizer, b: %f, sig_b: %d, log_b: %u\n", b, sl_b.sig_a, sl_b.log_a);
    double sum = inverse(add_sl(sl_a, sl_b));
    double sub = inverse(sub_sl(sl_a, sl_b));
    double mult = inverse(multi_sl(sl_a, sl_b));

    printf("SUM: %f, SUB: %f, MULTI: %f\n", sum, sub, mult);
}

void print_array_cplx_sl(const char * s, complex_sl_t a[], int n){
    printf("%s", s);
    for(int i=0; i<n; i++)  print_single_cplx_sl(a[i]);
    printf("\n");
}

void test_get_pow_e(){
    double ratio = 0.6;
    complex_sl_t t = get_pow_e(ratio);
    printf("MY: ");
    print_single_cplx_sl(t);
    
    printf("Default: ");
    cplx t2 = cexp(I * PI * ratio);
    if (!cimag(t2))
			printf("%g ", creal(t2));
		else
			printf("(%g, %g) ", creal(t2), cimag(t2));
}

/*
    Testing real inputs
*/
void test_log_fft(double* input_vec, double* output_vec, int n)
{
    complex_sl_t* buf = (complex_sl_t*) malloc(n * sizeof(complex_sl_t));

    for(int i=0; i<n; i++){
        buf[i].real = quantizer(input_vec[i]);
        buf[i].imag = quantizer(0);
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

void test_default_fft(double* input_vec, double* output_vec, int n)
{
    cplx* buf = (cplx*)malloc(n * sizeof(cplx));
    for(int i=0; i<n; i++){
        buf[i] = inverse(quantizer(input_vec[i]));
    }  
	//  cplx buf[8] = {0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0};
	
#ifdef _DEBUG_PRINT_
    show("Data: ", buf, n);
    printf("\n");
#endif

	fft(buf, n);

#ifdef _DEBUG_PRINT_
	show("FFT : ", buf, n);
    printf("\n");
#endif

    for(int i=0; i<n; i++)  output_vec[i] = creal(buf[i]);

    free(buf);
}

void test_n_point_fft(int n){
    //  Initialize random data
    double* input_vec = (double*)malloc(n * sizeof(double));
    double* default_output_vec = (double*)malloc(n * sizeof(double));
    double* sigNlog_output_vec = (double*)malloc(n * sizeof(double));
    double* variance_vec = (double*)malloc(n * sizeof(double));

    for(int i=0; i<n; i++)  input_vec[i] = randfrom(-1, 1);

    test_default_fft(input_vec, default_output_vec, n);
    test_log_fft(input_vec, sigNlog_output_vec, n);
    
    double error_sum = 0.0;
    double sum = 0.0;
    for(int i=0; i<n; i++){
        double noise = sigNlog_output_vec[i] - default_output_vec[i];
        double signal = default_output_vec[i];
        variance_vec[i] = (noise * noise);
        sum += signal * signal;
        error_sum += noise * noise;
    }

#ifdef _DEBUG_PRINT_
    printf("The variance: [");
    for(int i=0; i<n; i++)  printf("%g ", variance_vec[i]);
    printf("]\n");
#endif

    double average = error_sum / sum;

#ifdef _DEBUG_PRINT_
    printf("Average variance: %g, first element: %g, default: %g, new: %g\n", average, variance_vec[0], default_output_vec[0], sigNlog_output_vec[0]);
#endif

    printf("N: %d, SNR: %g\n", n, average);

    free(input_vec);
    free(default_output_vec);
    free(sigNlog_output_vec);
    free(variance_vec);
}

void simple_test(){
    double input_vec[8] = {0.250, -0.178, 0.511, 0.006, 0.557, -0.883, -0.374, -0.295};
    double output_vec1[8] = {1.0,1.0,1.0,1.0,0,0,0,0};
    double output_vec2[8] = {1.0,1.0,1.0,1.0,0,0,0,0};
    test_default_fft(input_vec, output_vec1, 8);
    test_log_fft(input_vec, output_vec2, 8);
}

int main(){
    //  simple_test();
    for(int i=1; i<12; i++){
        test_n_point_fft(1 << i);
    }
    
    
    return 0;
}
