#ifndef COMPLEX_MUL_H
#define COMPLEX_MUL_H

#include <complex>
#include <cmath>

using data_t = float;
using complex_t = std::complex<data_t>;

// Multiply two complex numbers (a_real + j*a_imag) * (b_real + j*b_imag).
void complex_mul(data_t a_real,
                 data_t a_imag,
                 data_t b_real,
                 data_t b_imag,
                 data_t &c_real,
                 data_t &c_imag);

#endif // COMPLEX_MUL_H
