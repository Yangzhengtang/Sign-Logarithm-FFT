#include <stdio.h>
#include <math.h>
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

int sg_T = 8;
int sg_B = (int)(sizeof(fixed_point_t))*8 - FIXED_POINT_FRACTIONAL_BITS;         //  Two parameters mentioned in the paper

typedef struct sig_log_t{
    char sig_a;
    fixed_point_t log_a;
}   sig_log_t;

sig_log_t multi_sl(sig_log_t a, sig_log_t b){
    fixed_point_t Kc = double_to_fixed(log2(sg_T));
    sig_log_t ret;
    ret.log_a = a.log_a + b.log_a - Kc;
    ret.sig_a = a.sig_a ^ b.sig_a;
    return ret;
}

/*
    Two helper functions used in the add and sub op.
    Possible Bug Here !!!
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
    if(a.log_a > b.log_a){
        ret.sig_a = a.sig_a;
        ret.log_a = a.log_a + gamma_func(b.log_a - a.log_a);
    }
    else{
        ret.sig_a = b.sig_a ? 0 : 1;
        ret.log_a = b.log_a + gamma_func(a.log_a - b.log_a);
    }
    return ret;
}

sig_log_t add_sl(sig_log_t a, sig_log_t b){
    if(a.sig_a == 0 && b.sig_a == 1){
        sig_log_t _b = b;
        _b.sig_a = 0;
        return _sub_sl(a, _b);
    }
    if(a.sig_a == 1 && b.sig_a == 0){
        sig_log_t _a = a;
        _a.sig_a = 0;
        return _sub_sl(b, _a);
    }
    return _add_sl(a, b);
}

sig_log_t sub_sl(sig_log_t a, sig_log_t b){
    if(a.sig_a == 0 && b.sig_a == 1){
        sig_log_t _b = b;
        _b.sig_a = 0;
        return _add_sl(a, _b);
    }
    if(a.sig_a == 1 && b.sig_a == 0){
        sig_log_t _a = a;
        _a.sig_a = 0;
        sig_log_t t = _add_sl(_a, b);
        t.sig_a = 1;
        return t;
    }
    if(a.sig_a == 1 && b.sig_a == 1){
        sig_log_t _a = a;
        _a.sig_a = 0;
        sig_log_t _b = b;
        _b.sig_a = 0;
        sig_log_t t = _add_sl(_a, _b);
        t.sig_a = 1;
        return t;
    }
    return _sub_sl(a, b);
}


typedef struct complex_sl_t{
    sig_log_t imag;
    sig_log_t real;
}   complex_sl_t;

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
        ret.sig_a = 1;
    else if(x > 0)
        ret.sig_a = 0;
    else{
        printf("Quantizing a 0: %f\n", x);  //  This case should be avoided.
        ret.sig_a = 1;
    }
    return ret;
}

double inverse(sig_log_t n){
    int sign = 1 - 2 * n.sig_a;
    double log_a_double = fixed_to_double(n.log_a);
    double value = (pow(2, log_a_double)) / sg_T;
    return (sign * value);
}

/*
    Xa = Xa + W * Xb
    Xb = Xa - W * Xa
    W = exp(−2πi/N k) Xk+N/2
*/
void log_butterfly_unit(complex_sl_t* xa, complex_sl_t* xb, complex_sl_t w){

}
