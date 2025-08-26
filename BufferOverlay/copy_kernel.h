#pragma once
#include <complex>
#ifndef N
#define N 64
#endif
static constexpr int TOTAL = N * N;
using complex_t = std::complex<float>;
extern "C" void copy_kernel(const complex_t* in, complex_t* out);
