// diff_losses.cpp
#include "diff_losses.h"

// ---- Implementación de wrappers (una firma, pragmas dentro) ----
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
  data_t inv = (data_t)1.0f / mag2;
  complex_t p = c_mul(c_conj(den), num);
  return complex_t(c_re(p)*inv, c_im(p)*inv);
}

static bool c_near_zero(const complex_t &z, data_t eps2) {
  #pragma HLS INLINE off
  return c_abs2(z) < eps2;
}

// --- Constantes numéricas (float) ---
static constexpr data_t Lx   = 45e-6f;
static constexpr data_t Ly   = 45e-6f;
static constexpr data_t dx   = Lx / N;
static constexpr data_t dy   = Ly / N;
static constexpr data_t dz   = 1e-4f;
static constexpr data_t k    = 7853981.6339f;   // ~2*pi/0.8e-6
static constexpr data_t eps  = 1e-12f;
static constexpr data_t eps2 = eps * eps;

static const complex_t C_ONE  = complex_t(1.0f, 0.0f);
static const complex_t C_ZERO = complex_t(0.0f, 0.0f);

// División compleja con protección (única firma)
static complex_t complex_div(const complex_t &num, const complex_t &den) {
  #pragma HLS INLINE off
  return c_div_safe(num, den, eps2);
}

// b = off*x0[i-1] + dp*x0[i] + off*x0[i+1], bordes dp1/dp2
static void compute_b_vector(
    complex_t dp,  complex_t dp1, complex_t dp2, complex_t off,
    const complex_t x0[N],
          complex_t b [N]
){
  #pragma HLS INLINE off
  #pragma HLS PIPELINE II=20   // evitar choques de puertos BRAM

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

// Thomas tridiagonal: (off, dp, off) con bordes dp1/dp2
static void thomas_solver(
    complex_t dp,  complex_t dp1, complex_t dp2, complex_t off,
    complex_t b[N],
    complex_t x[N]
){
  #pragma HLS INLINE off
  #pragma HLS PIPELINE II=20

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
  #pragma HLS PIPELINE II=20

  complex_t x0[N], b[N], x[N];

  complex_t ung   = c_make(0.0f, dz / (4 * k * dx * dx));
  complex_t neg2u = c_mul(c_make(-2.0f, 0.0f), ung);
  complex_t pos2u = c_mul(c_make( 2.0f, 0.0f), ung);

  for (int j = 0; j < N; ++j){
    #pragma HLS PIPELINE II=20

    // cargar columna j en x0
    for (int i = 0; i < N; ++i){
      #pragma HLS PIPELINE II=20
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
    thomas_solver(dp_A, dp1_A, dp2_A, c_scale(ung, -1.0f), b, x);

    // escribir salida
    for (int i = 0; i < N; ++i){
      #pragma HLS PIPELINE II=20
      out[i][j] = x[i];
    }
  }
}

// ADI en Y (filas i)
static void adi_y(const complex_t in[N][N], complex_t out[N][N]){
  #pragma HLS INLINE off
  #pragma HLS PIPELINE II=20

  complex_t x0[N], b[N], x[N];

  complex_t ung   = c_make(0.0f, dz / (4 * k * dy * dy));
  complex_t neg2u = c_mul(c_make(-2.0f, 0.0f), ung);
  complex_t pos2u = c_mul(c_make( 2.0f, 0.0f), ung);

  for (int i = 0; i < N; ++i){
    #pragma HLS PIPELINE II=20

    // cargar fila i en x0
    for (int j = 0; j < N; ++j){
      #pragma HLS PIPELINE II=20
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
    thomas_solver(dp_A, dp1_A, dp2_A, c_scale(ung, -1.0f), b, x);

    // escribir salida
    for (int j = 0; j < N; ++j){
      #pragma HLS PIPELINE II=20
      out[i][j] = x[j];
    }
  }
}

void diff_losses(
    const complex_t phi_in[N][N],
          complex_t phi_out[N][N],
    int steps
){
  #pragma HLS INTERFACE s_axilite port=return bundle=control
  #pragma HLS INTERFACE m_axi     port=phi_in  offset=slave bundle=gmem depth=4096
  #pragma HLS INTERFACE s_axilite port=phi_in  bundle=control
  #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=gmem depth=4096
  #pragma HLS INTERFACE s_axilite port=phi_out bundle=control
  #pragma HLS INTERFACE s_axilite port=steps  bundle=control

  // BRAM 2P, sin particiones (austero)
  complex_t phi[N][N];
  complex_t tmp[N][N];
  #pragma HLS BIND_STORAGE variable=phi type=ram_2p impl=uram
  #pragma HLS BIND_STORAGE variable=tmp type=ram_2p impl=uram

  // Copia de entrada (pipeline humilde)
  for (int i = 0; i < N; ++i){
    #pragma HLS PIPELINE II=8
    for (int j = 0; j < N; ++j){
      phi[i][j] = phi_in[i][j];
    }
  }

  // Núcleo ADI secuencial repetido steps veces
  for (int s = 0; s < steps; ++s){
    #pragma HLS PIPELINE II=20
    adi_x(phi, tmp);
    adi_y(tmp, phi);
  }

  // Copia de salida (pipeline humilde)
  for (int i = 0; i < N; ++i){
    #pragma HLS PIPELINE II=8
    for (int j = 0; j < N; ++j){
      phi_out[i][j] = phi[i][j];
    }
  }
}
