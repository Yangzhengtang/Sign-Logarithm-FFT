typedef short fixed_point_t;

typedef struct sig_log_t{
    int sig_a;
    fixed_point_t log_a;
}   sig_log_t;


typedef struct complex_sl_t{
    sig_log_t real;
    sig_log_t imag;
}   complex_sl_t;

sig_log_t quantizer(double x);
double inverse(sig_log_t n);
