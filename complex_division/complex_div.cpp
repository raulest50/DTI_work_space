#include "complex_div.h"

void complex_div(complex_t a, complex_t b, complex_t &c) {
#pragma HLS INTERFACE s_axilite port=a bundle=CTRL
#pragma HLS INTERFACE s_axilite port=b bundle=CTRL
#pragma HLS INTERFACE s_axilite port=c bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL

    c = a / b;
}
