#ifndef COMPLEX_DIV_H
#define COMPLEX_DIV_H

#include <complex>

using data_t = float;
using complex_t = std::complex<data_t>;

void complex_div(complex_t a, complex_t b, complex_t &c);

#endif // COMPLEX_DIV_H
