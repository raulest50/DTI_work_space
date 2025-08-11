// diffraction_only.h
#ifndef DIFFRACTION_ONLY_H
#define DIFFRACTION_ONLY_H

#include <complex>

const int N = 64;
using data_t    = float;
using complex_t = std::complex<data_t>;

void diffraction_only(const complex_t phi_in[N][N],
                      complex_t phi_out[N][N]);

#endif // DIFFRACTION_ONLY_H
