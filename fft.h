#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "log_fft.h"

#ifdef _SUPER_DEBUG_
    #define _DF_FFT_DEBUG_
#endif

double PI = 3.14159265358979323846;
typedef double complex cplx;
 
 
void show(const char * s, cplx buf[], int n) {
	printf("%s", s);
	for (int i = 0; i < n; i++)
		if (!cimag(buf[i]))
			printf("(%.3f, %.3f) ", creal(buf[i]), 0.0);
		else
			printf("(%.3f, %.3f) ", creal(buf[i]), cimag(buf[i]));
}

void _fft(cplx* buf, cplx* out, int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			
			
			complex_sl_t w = get_pow_e(-1.0 * i/n);
			cplx cplx_w = inverse(w.real) + I * inverse(w.imag);
			cplx t = cplx_w * out[i+step];
			
			//	cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
	#ifdef _DF_FFT_DEBUG_
        printf("_fft: n(%d), step(%d), the buffer: ", n, step);
        show("", buf, n);
        printf("\n");
    #endif
}
 
void fft(cplx* buf, int n)
{
	cplx* out = (cplx*) malloc(n * sizeof(cplx));
	for (int i = 0; i < n; i++) out[i] = buf[i];

	_fft(buf, out, n, 1);

	free(out);
}