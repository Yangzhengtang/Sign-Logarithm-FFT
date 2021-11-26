#include "log_fft.h"

void test_quantize_and_inverse(){
    double x = 0.314159263838438;
    char sig_x = 0;
    fixed_point_t  log_x = 0;
    sig_log_t sl_x = quantizer(x);

    printf("After quantizer, x: %f, sig_x: %d, log_x: %hu\n", x, sl_x.sig_a, sl_x.log_a);

    double new_x = inverse(sl_x);

    printf("new_x: %f\n", new_x);
}

void test_op(){
    double a = -0.75;
    double b = -0.25;
    sig_log_t sl_a = quantizer(a);
    sig_log_t sl_b = quantizer(b);

    printf("After quantizer, a: %f, sig_a: %d, log_a: %u\n", a, sl_a.sig_a, sl_a.log_a);
    printf("After quantizer, b: %f, sig_b: %d, log_b: %u\n", b, sl_b.sig_a, sl_b.log_a);
    double sum = inverse(add_sl(sl_a, sl_b));
    double sub = inverse(sub_sl(sl_a, sl_b));
    double mult = inverse(multi_sl(sl_a, sl_b));

    printf("SUM: %f, SUB: %f, MULTI: %f\n", sum, sub, mult);
}

int main(){
    test_op();
    return 0;
}