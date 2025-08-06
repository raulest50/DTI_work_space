#include "complex_mul.h"

void complex_mul(data_t a_real,
                 data_t a_imag,
                 data_t b_real,
                 data_t b_imag,
                 data_t &c_real,
                 data_t &c_imag) {
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE s_axilite port=a_real bundle=CTRL
#pragma HLS INTERFACE s_axilite port=a_imag bundle=CTRL
#pragma HLS INTERFACE s_axilite port=b_real bundle=CTRL
#pragma HLS INTERFACE s_axilite port=b_imag bundle=CTRL
#pragma HLS INTERFACE s_axilite port=c_real bundle=CTRL
#pragma HLS INTERFACE s_axilite port=c_imag bundle=CTRL

    complex_t a(a_real, a_imag);
    complex_t b(b_real, b_imag);
    complex_t c = a * b;
    c_real = c.real();
    c_imag = c.imag();
}
