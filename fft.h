#include <stdio.h>
#include <math.h>
#include <complex.h>
 
double PI = 3.14159265358979323846;
typedef double complex cplx;
 
void _fft(cplx* buf, cplx* out, int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
 
void fft(cplx* buf, int n)
{
	cplx* out = (cplx*) malloc(n * sizeof(cplx));
	for (int i = 0; i < n; i++) out[i] = buf[i];

	_fft(buf, out, n, 1);

	free(out);
}
 
 
void show(const char * s, cplx buf[], int n) {
	printf("%s", s);
	for (int i = 0; i < n; i++)
		if (!cimag(buf[i]))
			printf("(%.3f, %.3f) ", creal(buf[i]), 0.0);
		else
			printf("(%.3f, %.3f) ", creal(buf[i]), cimag(buf[i]));
}
 
 