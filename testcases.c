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