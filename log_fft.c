#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "utils.h"

/*
    Double precission number to sign/logarithm number
    input x should be in (0, 1)
*/
void quantizer(double x, int t, int b, char* sig_x, fixed_point_t* log_x){
    double log_x_double = 0;
    double trial = ((x * t) > 0) ? (x * t) : (x * t * -1);

    if(trial >= 1)
        log_x_double = log2(trial);   //  This should be a 17-bit logarithmic number, but don't know how many fraction bits it should contain.
    else
        log_x_double = 0;
    *log_x = double_to_fixed(log_x_double);

    if(x < 0)
        *sig_x = 1;
    else if(x > 0)
        *sig_x = 0;
    else{
        printf("Quantizing a 0: %f\n", x);  //  This case should be avoided.
        *sig_x = 1;
    }
        
}

void inverse(char sig_x, fixed_point_t log_x, int t, double* x){
    int sign = 1 - 2 * sig_x;
    double log_x_double = fixed_to_double(log_x);
    double value = (double)(pow(2, log_x_double)) / t;
    *x = sign * value;
}

void test_quantize_and_inverse(){
    double x = 0.314159263838438;
    char sig_x = 0;
    fixed_point_t  log_x = 0;
    quantizer(x, 8, 3, &sig_x, &log_x);

    printf("After quantizer, x: %f, sig_x: %d, log_x: %hu\n", x, sig_x, log_x);

    double new_x = 0.0;
    inverse(sig_x, log_x, 8, &new_x);

    printf("new_x: %f\n", new_x);
}





int main(){
    test_quantize_and_inverse();
    return 0;
}