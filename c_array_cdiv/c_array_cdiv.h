#ifndef C_ARRAY_CDIV_H
#define C_ARRAY_CDIV_H

#if __has_include(<hls_x_complex.h>)
#  include <hls_x_complex.h>
#  include <hls_stream.h>
#else
#  include "hls_stub.h"
#endif

const int N = 64;
using data_t    = float;
using complex_t = hls::x_complex<data_t>;

void c_array_cdiv(const complex_t in[N],
                  complex_t out[N],
                  data_t div_real,
                  data_t div_imag);

#endif // C_ARRAY_CDIV_H
