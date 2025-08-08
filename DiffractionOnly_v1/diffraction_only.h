// diffraction_only.h
#ifndef DIFFRACTION_ONLY_H
#define DIFFRACTION_ONLY_H

#if __has_include(<hls_x_complex.h>)
#  include <hls_x_complex.h>
#else
#  include "hls_stub.h"
#endif

const int N = 64;
using data_t    = float;
using complex_t = hls::x_complex<data_t>;

void diffraction_only(const complex_t phi_in[N][N],
                      complex_t phi_out[N][N]);

#endif // DIFFRACTION_ONLY_H
