#include "log_fft.h"
#include "fft.h"
//  utils is included from log_fft

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

void print_single_cplx_sl(complex_sl_t a){
    double imag = inverse(a.imag);
    double real = inverse(a.real);
    // double size = sqrt(imag * imag + real * real);

    printf("(%.3f, %.3f) ", real, imag);
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
int test_log_fft(double* input_vec, int n)
{
    complex_sl_t* buf = (complex_sl_t*) malloc(n * sizeof(complex_sl_t));

    for(int i=0; i<n; i++){
        buf[i].real = quantizer(input_vec[i]);
        buf[i].imag = quantizer(0);
    }

    print_array_cplx_sl("Data: ", buf, n);
    sl_fft(buf, n);
    print_array_cplx_sl("FFT : ", buf, n);
 
    free(buf);

	return 0;
}

int test_default_fft(double* input_vec, int n)
{
    cplx* buf = (cplx*)malloc(n * sizeof(cplx));
    for(int i=0; i<n; i++)  buf[i] = input_vec[i];
	//  cplx buf[8] = {0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0};
	show("Data: ", buf, n);
	fft(buf, 8);
	show("\nFFT : ", buf, n);
    printf("\n");

    free(buf);
	return 0;
}

void test_n_point_fft(int n){
    //  Initialize random data
    double* input_vec = (double*)malloc(n * sizeof(double));
    for(int i=0; i<n; i++)  input_vec[i] = randfrom(0, 1);

    test_default_fft(input_vec, n);
    test_log_fft(input_vec, n);
    
    free(input_vec);
}

void simple_test(){
    double input_vec[8] = {1.0,1.0,1.0,1.0,0,0,0,0};

    test_default_fft(input_vec, 8);
    test_log_fft(input_vec, 8);
}

int main(){
    test_n_point_fft(8);
    
    return 0;
}
