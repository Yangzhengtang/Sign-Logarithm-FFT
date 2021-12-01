#ifndef _LOGFFT_INCL_
#define _LOGFFT_INCL_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "log_fft.h"
#include <time.h>


#define FIXED_POINT_FRACTIONAL_BITS 11

double fixed_to_double(fixed_point_t input)
{
    return ((double)input / (double)(1 << FIXED_POINT_FRACTIONAL_BITS));
}

fixed_point_t double_to_fixed(double input)
{
    return (fixed_point_t)(round(input * (1 << FIXED_POINT_FRACTIONAL_BITS)));
}

static int random_initialized = 0;

/**
 * Generator of a random double value
 */
double randfrom(double min, double max)
{
    if(random_initialized == 0){
        srand(time(NULL));
        random_initialized = 1;
    }
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

#ifdef _SUPER_DEBUG_
    #define _LOG_FFT_DEBUG_
#endif

int sg_T = 1 << 5;
int sg_B = (int)(sizeof(fixed_point_t))*8 - FIXED_POINT_FRACTIONAL_BITS;         //  Two parameters mentioned in the paper

void copy_sl(sig_log_t dst, sig_log_t src){
    dst[0] = src[0];
    dst[1] = src[1];
}

void multi_sl(sig_log_t a, sig_log_t b, sig_log_t ret){
    fixed_point_t Kc = double_to_fixed(log2(sg_T));
    ret[0] = a[0] * b[0] - Kc;
    ret[1] = a[1] + b[1];
}

fixed_point_t beta_func(fixed_point_t x){
    double t = log2(pow(2, fixed_to_double(x)) + 1);
    return double_to_fixed(t);
}

fixed_point_t gamma_func(fixed_point_t x){
    double t = log2(1 - pow(2, fixed_to_double(x)));
    return double_to_fixed(t);
}

void _add_sl(sig_log_t a, sig_log_t b, sig_log_t ret){
    if(a[1] >= b[1]){
        ret[0] = a[0];
        ret[1] = a[1] + beta_func(b[1] - a[1]);
    }
    else{
        ret[0] = b[0];
        ret[1] = b[1] + beta_func(a[1] - b[1]);
    }
}

void _sub_sl(sig_log_t a, sig_log_t b, sig_log_t ret){
    if(a[1] == b[1]){     //  Don't know how to handle 0
        ret[0] = 0;
        ret[1] = 0;
        return;
    }
    if(a[1] > b[1]){
        ret[0] = a[0];
        ret[1] = a[1] + gamma_func(b[1] - a[1]);
    }
    else{
        ret[0] = - b[0];
        ret[1] = b[1] + gamma_func(a[1] - b[1]);
    }
    return;
}

void add_sl(sig_log_t a, sig_log_t b, sig_log_t ret){
    if(a[0] == 0){
        ret[0] = b[0];
        ret[1] = b[1];
        return;
    }
    if(b[0] == 0){
        ret[0] = a[0];
        ret[1] = a[1];
        return;
    }

    if(a[0] == 1 && b[0] == -1){
        sig_log_t _b;
        _b[0] = b[0];
        _b[1] = b[1];
        _b[0] = 1;
        _sub_sl(a, _b, ret);
        return;
    }
    if(a[0] == -1 && b[0] == 1){
        sig_log_t _a;
        _a[0] = a[0];
        _a[1] = a[1];
        _a[0] = 1;
        _sub_sl(b, _a, ret);
        return;
    }
    if(a[0] == -1 && b[0] == -1){
        sig_log_t _a;
        _a[0] = a[0];
        _a[1] = a[1];
        _a[0] = 1;

        sig_log_t _b;
        _b[0] = b[0];
        _b[1] = b[1];
        _b[0] = 1;
        _add_sl(_a, _b, ret);
        ret[0] = -1;
        return;
    }
    _add_sl(a, b, ret);
    return;
}

void sub_sl(sig_log_t a, sig_log_t b, sig_log_t ret){
	if(b[0] == 0){
    	ret[0] = a[0];
        ret[1] = a[1];
    	return;
    }
    if(a[0] == 0){
        ret[0] = b[0];
        ret[1] = b[1];
        ret[0] = -1 * (ret[0]);
        return;
    }
    
    if(a[0] == 1 && b[0] == -1){
        sig_log_t _b;
        _b[0] = b[0];
        _b[1] = b[1];
        _b[0] = 1;
        _add_sl(a, _b, ret);
        return;
    }
    if(a[0] == -1 && b[0] == 1){
        sig_log_t _a;
        _a[0] = a[0];
        _a[1] = a[1];
        _a[0] = 1;
        _add_sl(_a, b, ret);
        ret[0] = -1;
        return;
    }
    if(a[0] == -1 && b[0] == -1){
        sig_log_t _a;
        _a[0] = a[0];
        _a[1] = a[1];
        _a[0] = 1;

        sig_log_t _b;
        _b[0] = b[0];
        _b[1] = b[1];
        _b[0] = 1;
        _sub_sl(_b, _a, ret);
        return;
    }

	_sub_sl(a, b, ret);
    return;
}

/**
 * Complex sig_log number structure, and according operation:
 * 
 */

void copy_cplx_sl(complex_sl_t dst, complex_sl_t src){
    for(int i=0; i<4; i++)  dst[i] = src[i];
}

void copy_cplx_to_two_comp(sig_log_t real, sig_log_t imag, complex_sl_t src){
    real[0] = src[0];
    real[1] = src[1];
    imag[0] = src[2];
    imag[1] = src[3];
}

void copy_two_comp_to_cplx(complex_sl_t dst, sig_log_t real, sig_log_t imag){
    dst[0] = real[0];
    dst[1] = real[1];
    dst[2] = imag[0];
    dst[3] = imag[1];
}

void multi_cplx_sl(complex_sl_t a, complex_sl_t b, complex_sl_t ret){
    sig_log_t ret_real;
    sig_log_t ret_imag;

    sig_log_t a_real;
    sig_log_t a_imag;
    sig_log_t b_real;
    sig_log_t b_imag;

    copy_cplx_to_two_comp(a_real, a_imag, a);
    copy_cplx_to_two_comp(b_real, b_imag, b);

    sig_log_t temp_1;
    sig_log_t temp_2;

    //  ret.real = sub_sl(multi_sl(a.real, b.real), multi_sl(a.imag, b.imag));
    //  ret.imag = add_sl(multi_sl(a.real, b.imag), multi_sl(a.imag, b.real));

    multi_sl(a_real, b_real, temp_1);
    multi_sl(a_imag, b_imag, temp_2);
    sub_sl(temp_1, temp_2, ret_real);

    multi_sl(a_real, b_imag, temp_1);
    multi_sl(a_imag, b_real, temp_2);
    add_sl(temp_1, temp_2, ret_imag);

    copy_two_comp_to_cplx(ret, ret_real, ret_imag);
    return;
}

void add_cplx_sl(complex_sl_t a, complex_sl_t b, complex_sl_t ret){
    /*
    complex_sl_t ret;
    ret.real = add_sl(a.real, b.real);
    ret.imag = add_sl(a.imag, b.imag);
    */
    sig_log_t a_real;
    sig_log_t a_imag;
    sig_log_t b_real;
    sig_log_t b_imag;

    copy_cplx_to_two_comp(a_real, a_imag, a);
    copy_cplx_to_two_comp(b_real, b_imag, b);

    sig_log_t ret_real;
    sig_log_t ret_imag;

    add_sl(a_real, b_real, ret_real);
    add_sl(a_imag, b_imag, ret_imag);
    copy_two_comp_to_cplx(ret, ret_real, ret_imag);
    return;
}

void sub_cplx_sl(complex_sl_t a, complex_sl_t b, complex_sl_t ret){
    /*
    complex_sl_t ret;
    ret.real = sub_sl(a.real, b.real);
    ret.imag = sub_sl(a.imag, b.imag);
    return ret;
    */
    sig_log_t a_real;
    sig_log_t a_imag;
    sig_log_t b_real;
    sig_log_t b_imag;

    copy_cplx_to_two_comp(a_real, a_imag, a);
    copy_cplx_to_two_comp(b_real, b_imag, b);

    sig_log_t ret_real;
    sig_log_t ret_imag;

    sub_sl(a_real, b_real, ret_real);
    sub_sl(a_imag, b_imag, ret_imag);
    copy_two_comp_to_cplx(ret, ret_real, ret_imag);
    return;
}

/*
    Double precission number to sign/logarithm number
    input x should be in (0, 1)
*/
void quantizer(double x, sig_log_t ret){
    double log_a_double = 0;
    double trial = ((x * sg_T) > 0) ? (x * sg_T) : (x * sg_T * -1);

    if(trial >= 1)
        log_a_double = log2(trial);   //  This should be a 17-bit logarithmic number, but don't know how many fraction bits it should contain.
    else
        log_a_double = 0;
    ret[1] = double_to_fixed(log_a_double);

    if(x < 0)
        ret[0] = -1;
    else if(x > 0)
        ret[0] = 1;
    else{
        //  printf("Quantizing a 0: %f\n", x);  //  This case should be avoided.
        ret[0] = 0;
    }
    return;
}

double inverse(sig_log_t n){
    printf("on inverse\n");
    int sign = n[0];
    double log_a_double = fixed_to_double(n[1]);
    double value = (pow(2, log_a_double)) / sg_T;
    return (sign * value);
}

void print_single_cplx_sl(complex_sl_t a){

    sig_log_t a_real;
    sig_log_t a_imag;
    copy_cplx_to_two_comp(a_real, a_imag, a);

    double imag = inverse(a_imag);
    double real = inverse(a_real);

    printf("(%.3f, %.3f) ", real, imag);
}

void get_pow_e(double ratio, complex_sl_t ret){
    //  const double _PI = atan2(1, 1) * 4;
    const double _PI = 3.14159265358979323846;
    sig_log_t ret_real;
    sig_log_t ret_imag;
    quantizer(cos(_PI * ratio), ret_real);
    quantizer(sin(_PI * ratio), ret_imag);
    copy_two_comp_to_cplx(ret, ret_real, ret_imag);
    return;
}

#define FFT_POINT 8
void _sl_fft(complex_sl_t buf[FFT_POINT], complex_sl_t out[FFT_POINT], int step)
{
    printf("On %d\n", step);
	if (step < FFT_POINT) {
		_sl_fft(out, buf, step * 2);
		_sl_fft(out + step * 4, buf + step * 4, step * 2);
 
		for (int i = 0; i < FFT_POINT; i += 2 * step) {
            complex_sl_t w;
            get_pow_e(-1.0 * i / FFT_POINT, w);
            complex_sl_t temp;
            complex_sl_t temp_operand;
            //copy_cplx_sl(temp_operand, out[])
            multi_cplx_sl(w, out[(i+step)], temp);
			add_cplx_sl(out[i], temp, buf[(i/2)]);
			sub_cplx_sl(out[i], temp, buf[((i + FFT_POINT)/2)]);
		}
	}
}

/*
void _sl_fft_16(complex_sl_t* buf, complex_sl_t* out)
{
    return;
}

void _sl_fft_8(complex_sl_t* buf, complex_sl_t* out)
{
  if (8 < n) {
		_sl_fft_16(out, buf, n);
		_sl_fft_16(out + 8, buf + 8, n);

		for (int i = 0; i < n; i += 8) {
            complex_sl_t w = get_pow_e(-1.0 * i / n);
            complex_sl_t temp = multi_cplx_sl(w, out[i+8]);
			buf[i / 2]     = add_cplx_sl(out[i], temp);
			buf[(i + n)/2] = sub_cplx_sl(out[i], temp);
		}
	}
}

void _sl_fft_4(complex_sl_t* buf, complex_sl_t* out)
{
  if (4 < n) {
		_sl_fft_8(out, buf, n);
		_sl_fft_8(out + 4, buf + 4, n);

		for (int i = 0; i < n; i += 8) {
            complex_sl_t w = get_pow_e(-1.0 * i / n);
            complex_sl_t temp = multi_cplx_sl(w, out[i+4]);
			buf[i / 2]     = add_cplx_sl(out[i], temp);
			buf[(i + n)/2] = sub_cplx_sl(out[i], temp);
		}
	}
}


void _sl_fft_2(complex_sl_t* buf, complex_sl_t* out)
{
  if (2 < n) {
		_sl_fft_4(out, buf, n);
		_sl_fft_4(out + 2, buf + 2, n);

		for (int i = 0; i < n; i += 4) {
            complex_sl_t w = get_pow_e(-1.0 * i / n);
            complex_sl_t temp = multi_cplx_sl(w, out[i+2]);
			buf[i / 2]     = add_cplx_sl(out[i], temp);
			buf[(i + n)/2] = sub_cplx_sl(out[i], temp);
		}
	}
}



//  step = 1
void _sl_fft_1(complex_sl_t buf[FFT_POINT], complex_sl_t out[FFT_POINT])
{
  if (1 < FFT_POINT) {
		_sl_fft_2(out, buf);
		_sl_fft_2(out + 1 * 4, buf + 1 * 4);

		for (int i = 0; i < FFT_POINT; i += 2) {
            complex_sl_t w;
            get_pow_e(-1.0 * i / FFT_POINT, w);
            complex_sl_t temp;
            multi_cplx_sl(w, out[i+1], temp);
			add_cplx_sl(out[i * 4], temp, buf[(i / 2) * 4]);
			sub_cplx_sl(out[i * 4], temp, buf[((i + FFT_POINT)/2) * 4]);
		}
	}
}*/
 
void sl_fft(complex_sl_t buf[FFT_POINT])
{
	//complex_sl_t* out = (complex_sl_t*)malloc(n * sizeof(complex_sl_t));

	complex_sl_t out[FFT_POINT];

	for (int i = 0; i < FFT_POINT; i++)
        copy_cplx_sl(out[i], buf[i]);
 
    _sl_fft(buf, out, 1);
	//  _sl_fft_1(buf, out);
}

#endif
