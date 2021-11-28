#ifndef _LOGFFT_INCL_
#define _LOGFFT_INCL_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"

/*
    X0,...,N−1 ← ditfft2(x, N, s):             DFT of (x0, xs, x2s, ..., x(N-1)s):
    if N = 1 then
        X0 ← x0                                      trivial size-1 DFT base case
    else
        X0,...,N/2−1 ← ditfft2(x, N/2, 2s)             DFT of (x0, x2s, x4s, ..., x(N-2)s)
        XN/2,...,N−1 ← ditfft2(x+s, N/2, 2s)           DFT of (xs, xs+2s, xs+4s, ..., x(N-1)s)
        for k = 0 to N/2−1 do                        combine DFTs of two halves into full DFT:
            p ← Xk
            q ← exp(−2πi/N k) Xk+N/2
            Xk ← p + q 
            Xk+N/2 ← p − q
        end for
    end if
*/
//  #define _SUPER_DEBUG_

#ifdef _SUPER_DEBUG_
    #define _LOG_FFT_DEBUG_
#endif

int sg_T = 0b1000000000;
int sg_B = (int)(sizeof(fixed_point_t))*8 - FIXED_POINT_FRACTIONAL_BITS;         //  Two parameters mentioned in the paper

typedef struct sig_log_t{
    int sig_a;
    fixed_point_t log_a;
}   sig_log_t;

sig_log_t multi_sl(sig_log_t a, sig_log_t b){
    fixed_point_t Kc = double_to_fixed(log2(sg_T));
    sig_log_t ret;
    ret.log_a = a.log_a + b.log_a - Kc;
    ret.sig_a = a.sig_a * b.sig_a;
    return ret;
}

/*
    Two helper functions used in the add and sub op.
    Possible Bug Here !!!
    Now they only accept positive input.
*/
fixed_point_t beta_func(fixed_point_t x){
    double t = log2(pow(2, fixed_to_double(x)) + 1);
    return double_to_fixed(t);
}

fixed_point_t gamma_func(fixed_point_t x){
    double t = log2(1 - pow(2, fixed_to_double(x)));
    return double_to_fixed(t);
}

sig_log_t _add_sl(sig_log_t a, sig_log_t b){
    sig_log_t ret;
    if(a.log_a >= b.log_a){
        ret.sig_a = a.sig_a;
        ret.log_a = a.log_a + beta_func(b.log_a - a.log_a);
    }
    else{
        ret.sig_a = b.sig_a;
        ret.log_a = b.log_a + beta_func(a.log_a - b.log_a);
    }
    return ret;
}

sig_log_t _sub_sl(sig_log_t a, sig_log_t b){
    sig_log_t ret;
    if(a.log_a == b.log_a){     //  Don't know how to handle 0
        ret.sig_a = 0;
        ret.log_a = 0;
        return ret;
    }
    if(a.log_a > b.log_a){
        ret.sig_a = a.sig_a;
        ret.log_a = a.log_a + gamma_func(b.log_a - a.log_a);
    }
    else{
        ret.sig_a = - b.sig_a;
        ret.log_a = b.log_a + gamma_func(a.log_a - b.log_a);
    }
    return ret;
}

sig_log_t add_sl(sig_log_t a, sig_log_t b){
    if(a.sig_a == 0)    return b;
    if(b.sig_a == 0)    return a;
    if(a.sig_a == 1 && b.sig_a == -1){
        sig_log_t _b = b;
        _b.sig_a = 1;
        return _sub_sl(a, _b);
    }
    if(a.sig_a == -1 && b.sig_a == 1){
        sig_log_t _a = a;
        _a.sig_a = 1;
        return _sub_sl(b, _a);
    }
    if(a.sig_a == -1 && b.sig_a == -1){
        sig_log_t _a = a;
        _a.sig_a = 1;
        sig_log_t _b = b;
        _b.sig_a = 1;
        sig_log_t t = _add_sl(_a, _b);
        t.sig_a = -1;
        return t;
    }
    return _add_sl(a, b);
}

sig_log_t sub_sl(sig_log_t a, sig_log_t b){
    if(b.sig_a == 0)    return a;
    if(a.sig_a == 0){
        sig_log_t _b = b;
        _b.sig_a = -1 * (_b.sig_a);
        return _b;
    }
    
    if(a.sig_a == 1 && b.sig_a == -1){
        sig_log_t _b = b;
        _b.sig_a = 1;
        return _add_sl(a, _b);
    }
    if(a.sig_a == -1 && b.sig_a == 1){
        sig_log_t _a = a;
        _a.sig_a = 1;
        sig_log_t t = _add_sl(_a, b);
        t.sig_a = -1;
        return t;
    }
    if(a.sig_a == -1 && b.sig_a == -1){
        sig_log_t _a = a;
        _a.sig_a = 1;
        sig_log_t _b = b;
        _b.sig_a = 1;
        sig_log_t t = _sub_sl(_b, _a);
        return t;
    }
    return _sub_sl(a, b);
}

/**
 * Complex sig_log number structure, and according operation:
 * 
 */
typedef struct complex_sl_t{
    sig_log_t real;
    sig_log_t imag;
}   complex_sl_t;

complex_sl_t multi_cplx_sl(complex_sl_t a, complex_sl_t b){
    complex_sl_t ret;
    ret.real = sub_sl(multi_sl(a.real, b.real), multi_sl(a.imag, b.imag));
    ret.imag = add_sl(multi_sl(a.real, b.imag), multi_sl(a.imag, b.real));
    return ret;
}

complex_sl_t add_cplx_sl(complex_sl_t a, complex_sl_t b){
    complex_sl_t ret;
    ret.real = add_sl(a.real, b.real);
    ret.imag = add_sl(a.imag, b.imag);
    return ret;
}

complex_sl_t sub_cplx_sl(complex_sl_t a, complex_sl_t b){
    complex_sl_t ret;
    ret.real = sub_sl(a.real, b.real);
    ret.imag = sub_sl(a.imag, b.imag);
    return ret;
}

/*
    Double precission number to sign/logarithm number
    input x should be in (0, 1)
*/
sig_log_t quantizer(double x){
    sig_log_t ret;
    double log_a_double = 0;
    double trial = ((x * sg_T) > 0) ? (x * sg_T) : (x * sg_T * -1);

    if(trial >= 1)
        log_a_double = log2(trial);   //  This should be a 17-bit logarithmic number, but don't know how many fraction bits it should contain.
    else
        log_a_double = 0;
    ret.log_a = double_to_fixed(log_a_double);

    if(x < 0)
        ret.sig_a = -1;
    else if(x > 0)
        ret.sig_a = 1;
    else{
        //  printf("Quantizing a 0: %f\n", x);  //  This case should be avoided.
        ret.sig_a = 0;
    }
    return ret;
}

double inverse(sig_log_t n){
    int sign = n.sig_a;
    double log_a_double = fixed_to_double(n.log_a);
    double value = (pow(2, log_a_double)) / sg_T;
    return (sign * value);
}

void print_single_cplx_sl(complex_sl_t a){
    double imag = inverse(a.imag);
    double real = inverse(a.real);

    printf("(%.3f, %.3f) ", real, imag);
}

/*
    Xa = Xa + W * Xb
    Xb = Xa - W * Xb
    W = exp(−2πi/N k) Xk+N/2
*/
void log_butterfly_unit(complex_sl_t xa, complex_sl_t xb, complex_sl_t w){
    complex_sl_t temp = multi_cplx_sl(w, xb);
    xa = add_cplx_sl(xa, temp);
    xb = sub_cplx_sl(xa, temp);
}

complex_sl_t get_pow_e(double ratio){
    //  const double _PI = atan2(1, 1) * 4;
    const double _PI = 3.14159265358979323846;
    complex_sl_t ret;
    ret.real = quantizer(cos(_PI * ratio));
    ret.imag = quantizer(sin(_PI * ratio));
    return ret;
}


void _sl_fft(complex_sl_t* buf, complex_sl_t* out, int n, int step)
{
	if (step < n) {
		_sl_fft(out, buf, n, step * 2);
		_sl_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
            complex_sl_t w = get_pow_e(-1.0 * i / n);
            complex_sl_t temp = multi_cplx_sl(w, out[i+step]);
			buf[i / 2]     = add_cplx_sl(out[i], temp);
			buf[(i + n)/2] = sub_cplx_sl(out[i], temp);
#ifdef _LOG_FFT_DEBUG_
            if(n == 8 && step == 4 && i==0){
                printf("*** At debug point: buf[4]: ");
                print_single_cplx_sl(buf[(i + n)/2]);
                printf(", out[i]: ");
                print_single_cplx_sl(out[i]);
                printf(", w: ");
                print_single_cplx_sl(w);
                printf(", out[i+step]: ");
                print_single_cplx_sl(out[i+step]);
                printf(", temp: ");
                print_single_cplx_sl(temp);
                printf("***\n");
            }
#endif
		}
	}

    #ifdef _LOG_FFT_DEBUG_
        printf("_sl_fft: n(%d), step(%d), the buffer: ", n, step);
        for(int i=0; i<n; i++){
            double imag = inverse(buf[i].imag);
            double real = inverse(buf[i].real);
            printf("(%.3f, %.3f) ", real, imag);
        }
        printf("\n");
    #endif
}
 
void sl_fft(complex_sl_t* buf, int n)
{
	complex_sl_t* out = (complex_sl_t*)malloc(n * sizeof(complex_sl_t));

	for (int i = 0; i < n; i++)
        out[i] = buf[i];
 
	_sl_fft(buf, out, n, 1);

    free(out);
}

#endif