#pragma once
#include <complex>
#include "ap_int.h"
#include "ap_fixed.h"
#include "hls_stream.h"

#ifndef N
#define N 64
#endif

// Tipo base y complejo
using data_t    = float;
using complex_t = std::complex<float>;

// ---- Wrappers SIN plantillas (una única firma) ----
// Importante: evitar llamadas directas a operadores de std::complex
// en múltiples contextos. Centralizamos aquí y desactivamos inline.

static inline complex_t c_make(data_t re=0.f, data_t im=0.f) {
    return complex_t(re, im);
}
static inline data_t    c_re(const complex_t &z) { return z.real(); }
static inline data_t    c_im(const complex_t &z) { return z.imag(); }

static complex_t c_add (const complex_t &a, const complex_t &b);
static complex_t c_sub (const complex_t &a, const complex_t &b);
static complex_t c_scale(const complex_t &a, data_t s);
static complex_t c_mul (const complex_t &a, const complex_t &b);
static data_t    c_abs2(const complex_t &a);
static complex_t c_conj(const complex_t &a);
static complex_t c_div_safe(const complex_t &num, const complex_t &den, data_t eps2);
static bool      c_near_zero(const complex_t &z, data_t eps2);

// Top
void diff_losses(
    const complex_t phi_in[N][N],
          complex_t phi_out[N][N],
    int steps,
    data_t alpha,
    data_t beta
);
