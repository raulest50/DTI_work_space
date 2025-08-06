#ifndef THOMAS_SOLVER_H
#define THOMAS_SOLVER_H

#if __has_include(<hls_x_complex.h>)
#  include <hls_x_complex.h>
#  include <hls_stream.h>
#else
#  include "hls_stub.h"
#endif

#include <cmath>

const int N = 64;
using data_t = float;
using complex_t = hls::x_complex<data_t>;

// Solve a tridiagonal system with constant diagonals.  All coefficients and
// the right hand side are complex floating point values using single
// precision for the real and imaginary parts.  The main diagonal is
// ``dp`` everywhere except for the first and last elements that are ``dp1``
// and ``dp2`` respectively.  ``off`` is used for both sub and super
// diagonals.  ``b`` contains the right hand side vector and ``x`` will hold
// the solution.
void thomas_solver(complex_t dp,
                   complex_t dp1,
                   complex_t dp2,
                   complex_t off,
                   complex_t b[N],
                   complex_t x[N]);

#endif // THOMAS_SOLVER_H
