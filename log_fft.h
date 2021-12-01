#define FFT_POINT 8

/* typedef struct sig_log_t{
    int sig_a;
    fixed_point_t log_a;
}   sig_log_t; */

/* 
typedef struct complex_sl_t{
    sig_log_t real;
    sig_log_t imag;
}   complex_sl_t;
 */

typedef short fixed_point_t;

//  [0] : sig_a
//  [1] : log_a
typedef fixed_point_t sig_log_t[2];

//  [0] : real.sig_a
//  [1] : real.log_a
//  [2] : imag.sig_a
//  [3] : imag.log_a

typedef fixed_point_t complex_sl_t[4];

void quantizer(double x, sig_log_t ret);
double inverse(sig_log_t n);
