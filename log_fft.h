typedef short fixed_point_t;

/* typedef struct sig_log_t{
    int sig_a;
    fixed_point_t log_a;
}   sig_log_t; */


typedef fixed_point_t sig_log_t[2];
/* 
typedef struct complex_sl_t{
    sig_log_t real;
    sig_log_t imag;
}   complex_sl_t;
 */
typedef complex_sl_t fixed_point_t[2];

void quantizer(double x, sig_log_t ret);
double inverse(sig_log_t n);
