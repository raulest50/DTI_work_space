// diff_losses.cpp
#include "diff_losses.h"
#include <cmath>
#include <cstdio>  // For printf during csim (ignored in synthesis)

// Wrapper for exp() that works with ap_fixed types
// Converts to float for the exp calculation (compatible with both csim and synthesis)
static data_t fp_exp(data_t x) {
    #pragma HLS INLINE off
    // Convert to float, compute exp, convert back
    // expf() resolves to std::expf for csim, hls::expf for synthesis
    float x_f = (float)x;
    float result_f = expf(x_f);
    return (data_t)result_f;
}

// ---- Wrappers con firma única (evita clonación de std::complex) ----
static complex_t c_add(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)+c_re(b), c_im(a)+c_im(b));
}
static complex_t c_sub(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)-c_re(b), c_im(a)-c_im(b));
}
static complex_t c_scale(const complex_t &a, data_t s) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)*s, c_im(a)*s);
}
static complex_t c_mul(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  data_t ar = c_re(a), ai = c_im(a);
  data_t br = c_re(b), bi = c_im(b);
  return complex_t(ar*br - ai*bi, ar*bi + ai*br);
}
static data_t c_abs2(const complex_t &a) {
  #pragma HLS INLINE off
  data_t ar = c_re(a), ai = c_im(a);
  return ar*ar + ai*ai;
}
static complex_t c_conj(const complex_t &a) {
  #pragma HLS INLINE off
  return complex_t(c_re(a), -c_im(a));
}
static complex_t c_div_safe(const complex_t &num, const complex_t &den, data_t eps2) {
  #pragma HLS INLINE off
  data_t mag2 = c_abs2(den);
  if (mag2 < eps2) mag2 = eps2;
  data_t inv = data_t(1) / mag2;
  complex_t p = c_mul(c_conj(den), num);
  return complex_t(c_re(p)*inv, c_im(p)*inv);
}
static bool c_near_zero(const complex_t &z, data_t eps2) {
  #pragma HLS INLINE off
  return c_abs2(z) < eps2;
}

// --- Constantes numéricas ---
// Note: k=7853981.6339 is too large for ap_fixed<27,8>, so we use precomputed coefficients
// Physical parameters (kept as comments for reference):
//   Lx = Ly = 45e-6 m, dx = dy = Lx/N, dz = 1e-4 m, k = 2*pi/0.8e-6 ~ 7.85e6

// Precomputed coefficient: dz / (4 * k * dx * dx) where dx = 45e-6/64
// = 1e-4 / (4 * 7853981.6339 * (7.03125e-7)^2) = 6.4368
static const data_t UNG_COEFF = data_t(6.4368);

// dz / 4 for absorption calculations = 1e-4 / 4 = 2.5e-5
static const data_t DZ_OVER_4 = data_t(2.5e-5);

// Small epsilon for numerical stability
// IMPORTANT: ap_fixed<27,8> resolution is 2^-19 ≈ 1.9e-6
// eps values must be larger than this to be representable!
static const data_t eps  = data_t(1e-3);   // For general comparisons
static const data_t eps2 = data_t(1e-5);   // For division protection (must be > 2^-19)

static const complex_t C_ONE  = complex_t(data_t(1), data_t(0));
static const complex_t C_ZERO = complex_t(data_t(0), data_t(0));

// Absorption operators -------------------------------------------------
static void half_linear_absorption(complex_t phi[N][N], data_t alpha) {
  #pragma HLS INLINE off
  // exp(-alpha * dz/4) using precomputed DZ_OVER_4
  data_t arg = -alpha * DZ_OVER_4;
  data_t c = fp_exp(arg);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      phi[i][j] = c_scale(phi[i][j], c);
}

static void half_2photon_absorption(complex_t phi[N][N], data_t beta) {
  #pragma HLS INLINE off
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      // exp(-beta * dz/4 * |phi|^2) using precomputed DZ_OVER_4
      data_t intensity = c_abs2(phi[i][j]);
      data_t arg = -beta * DZ_OVER_4 * intensity;
      data_t atten = fp_exp(arg);
      phi[i][j] = c_scale(phi[i][j], atten);
    }
}

static void apply_N_midpoint(complex_t phi[N][N], data_t alpha, data_t beta) {
  #pragma HLS INLINE off
  half_linear_absorption(phi, alpha);
  half_2photon_absorption(phi, beta);
  half_2photon_absorption(phi, beta);
  half_linear_absorption(phi, alpha);
}

// División compleja protegida (una sola implementación)
static complex_t complex_div(const complex_t &num, const complex_t &den) {
  #pragma HLS INLINE off
  return c_div_safe(num, den, eps2);
}

// b = off*x0[i-1] + dp*x0[i] + off*x0[i+1], con dp1/dp2 en bordes
static void compute_b_vector(
    complex_t dp,  complex_t dp1, complex_t dp2, complex_t off,
    const complex_t x0[N],
          complex_t b [N]
){
  #pragma HLS INLINE off
  // Sin PIPELINE: evita presión de puertos y evita inlining/replicación
  b[0] = c_add( c_mul(x0[0], dp1), c_mul(x0[1], off) );
  for (int i = 1; i < N-1; ++i) {
    complex_t t = C_ZERO;
    t = c_add(t, c_mul(x0[i-1], off));
    t = c_add(t, c_mul(x0[i],   dp));
    t = c_add(t, c_mul(x0[i+1], off));
    b[i] = t;
  }
  b[N-1] = c_add( c_mul(x0[N-2], off), c_mul(x0[N-1], dp2) );
}

// Thomas tridiagonal: (off, dp, off) y bordes dp1/dp2
static void thomas_solver(
    complex_t dp,  complex_t dp1, complex_t dp2, complex_t off,
    complex_t b[N],
    complex_t x[N]
){
  #pragma HLS INLINE off

  complex_t c_prime[N];
  complex_t d_prime[N];

  // i = 0
  complex_t inv1 = complex_div(C_ONE, dp1);
  c_prime[0]     = c_mul(inv1, off);
  d_prime[0]     = c_mul(b[0], inv1);

  // forward
  for (int i = 1; i < N-1; ++i){
    complex_t denom = c_sub(dp, c_mul(c_prime[i-1], off));
    complex_t invd  = complex_div(C_ONE, denom);
    c_prime[i]      = c_mul(invd, off);
    d_prime[i]      = c_mul( c_sub(b[i], c_mul(d_prime[i-1], off)), invd );
  }

  // último
  {
    complex_t denom   = c_sub(dp2, c_mul(c_prime[N-2], off));
    complex_t invd    = complex_div(C_ONE, denom);
    d_prime[N-1]      = c_mul( c_sub(b[N-1], c_mul(d_prime[N-2], off)), invd );
  }

  // backward
  x[N-1] = d_prime[N-1];
  for (int i = N-2; i >= 0; --i){
    x[i] = c_sub(d_prime[i], c_mul(c_prime[i], x[i+1]));
  }
}

// ADI en X (columnas j)
static void adi_x(const complex_t in[N][N], complex_t out[N][N]){
  #pragma HLS INLINE off

  complex_t x0[N], b[N], x[N];

  // ung = i * dz/(4*k*dx^2), using precomputed UNG_COEFF
  complex_t ung   = c_make(data_t(0), UNG_COEFF);
  complex_t neg2u = c_mul(c_make(data_t(-2), data_t(0)), ung);
  complex_t pos2u = c_mul(c_make(data_t(2), data_t(0)), ung);

  for (int j = 0; j < N; ++j){
    // cargar columna j en x0
    for (int i = 0; i < N; ++i){
      x0[i] = in[i][j];
    }

    complex_t ratio0 = c_near_zero(x0[1],   eps2) ? C_ONE : complex_div(x0[0],   x0[1]);
    complex_t ratioN = c_near_zero(x0[N-2], eps2) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

    // Paso B
    complex_t dp_B  = c_add(neg2u, C_ONE);
    complex_t dp1_B = c_add(dp_B, c_mul(ung, ratio0));
    complex_t dp2_B = c_add(dp_B, c_mul(ung, ratioN));
    compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

    // Paso A
    complex_t dp_A  = c_add(pos2u, C_ONE);
    complex_t dp1_A = c_sub(dp_A, c_mul(ung, ratio0));
    complex_t dp2_A = c_sub(dp_A, c_mul(ung, ratioN));
    thomas_solver(dp_A, dp1_A, dp2_A, c_scale(ung, data_t(-1)), b, x);

    // escribir salida
    for (int i = 0; i < N; ++i){
      out[i][j] = x[i];
    }
  }
}

// ADI en Y (filas i)
static void adi_y(const complex_t in[N][N], complex_t out[N][N]){
  #pragma HLS INLINE off

  complex_t x0[N], b[N], x[N];

  // ung = i * dz/(4*k*dy^2), using precomputed UNG_COEFF (same as X since dx=dy)
  complex_t ung   = c_make(data_t(0), UNG_COEFF);
  complex_t neg2u = c_mul(c_make(data_t(-2), data_t(0)), ung);
  complex_t pos2u = c_mul(c_make(data_t(2), data_t(0)), ung);

  for (int i = 0; i < N; ++i){
    // cargar fila i en x0
    for (int j = 0; j < N; ++j){
      x0[j] = in[i][j];
    }

    complex_t ratio0 = c_near_zero(x0[1],   eps2) ? C_ONE : complex_div(x0[0],   x0[1]);
    complex_t ratioN = c_near_zero(x0[N-2], eps2) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

    // Paso B
    complex_t dp_B  = c_add(neg2u, C_ONE);
    complex_t dp1_B = c_add(dp_B, c_mul(ung, ratio0));
    complex_t dp2_B = c_add(dp_B, c_mul(ung, ratioN));
    compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

    // Paso A
    complex_t dp_A  = c_add(pos2u, C_ONE);
    complex_t dp1_A = c_sub(dp_A, c_mul(ung, ratio0));
    complex_t dp2_A = c_sub(dp_A, c_mul(ung, ratioN));
    thomas_solver(dp_A, dp1_A, dp2_A, c_scale(ung, data_t(-1)), b, x);

    // escribir salida
    for (int j = 0; j < N; ++j){
      out[i][j] = x[j];
    }
  }
}

void diff_losses(
    const complex_t phi_in[N][N],
          complex_t phi_out[N][N],
    int steps,
    data_t alpha,
    data_t beta
){
  #pragma HLS INTERFACE s_axilite port=return bundle=control
  #pragma HLS INTERFACE m_axi     port=phi_in  offset=slave bundle=gmem depth=4096
  #pragma HLS INTERFACE s_axilite port=phi_in  bundle=control
  #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=gmem depth=4096
  #pragma HLS INTERFACE s_axilite port=phi_out bundle=control
  #pragma HLS INTERFACE s_axilite port=steps  bundle=control
  #pragma HLS INTERFACE s_axilite port=alpha  bundle=control
  #pragma HLS INTERFACE s_axilite port=beta   bundle=control

  // Matrices internas en BRAM 2P (pequeño y rápido de sintetizar)
  complex_t phi[N][N];
  complex_t tmp[N][N];
  #pragma HLS BIND_STORAGE variable=phi type=ram_2p impl=bram
  #pragma HLS BIND_STORAGE variable=tmp type=ram_2p impl=bram

  // Copia de entrada (sin pipeline)
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      phi[i][j] = phi_in[i][j];

  // Repite ADI 'steps' veces (sin pipeline alrededor)
  for (int s = 0; s < steps; ++s){
    // Progress indicator for csim (ignored in synthesis)
    if (s % 20 == 0) {
      printf("[csim] Step %d / %d\n", s, steps);
      fflush(stdout);
    }
    adi_x(phi, tmp);
    apply_N_midpoint(tmp, alpha, beta);
    adi_y(tmp, phi);
  }
  printf("[csim] Completed all %d steps\n", steps);

  // Copia de salida (sin pipeline)
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      phi_out[i][j] = phi[i][j];
}
