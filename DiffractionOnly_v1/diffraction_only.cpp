// diffraction_only.cpp
#include "diffraction_only.h"

// keep constants simple; constexpr is fine and free
static constexpr data_t Lx   = 45e-6f;
static constexpr data_t Ly   = 45e-6f;
static constexpr data_t dx   = Lx / N;
static constexpr data_t dy   = Ly / N;
static constexpr data_t dz   = 1e-4f;
static constexpr data_t k    = 7853981.6339f;   // ~2*pi/0.8e-6
static constexpr data_t eps  = 1e-12f;
static constexpr data_t eps2 = eps * eps;
static const     complex_t C_ONE{1.0f, 0.0f};

inline bool is_near_zero(const complex_t &z) {
    data_t re = z.real();
    data_t im = z.imag();
    return (re*re + im*im) < eps2;
}

// Minimal manual complex divide (no operator/), with tiny guard
#pragma HLS INLINE
static complex_t complex_div(const complex_t &num, const complex_t &den) {
    data_t re = den.real();
    data_t im = den.imag();
    data_t mag2 = re*re + im*im;
    if (mag2 < eps2) mag2 = eps2;
    data_t inv = (data_t)1.0f / mag2;
    complex_t cd(re, -im);
    complex_t p  = cd * num;
    return complex_t(p.real() * inv, p.imag() * inv);
}

static void compute_b_vector(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t off,
    const complex_t x0[N],
          complex_t b[N]
) {
    b[0] = x0[0] * dp1 + x0[1] * off;
    for (int i = 1; i < N-1; i++) {
        b[i] = x0[i-1] * off + x0[i] * dp + x0[i+1] * off;
    }
    b[N-1] = x0[N-2] * off + x0[N-1] * dp2;
}

static void thomas_solver(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t off,
          complex_t b[N],
          complex_t x[N]
) {
    complex_t c_prime[N];
    complex_t d_prime[N];

    // first iteration
    complex_t inv1 = complex_div(C_ONE, dp1);
    c_prime[0]     = inv1 * off;
    d_prime[0]     = b[0] * inv1;

    // forward sweep
    for (int i = 1; i < N-1; i++) {
        complex_t denom = dp - c_prime[i-1] * off;
        complex_t invd  = complex_div(C_ONE, denom);
        c_prime[i]      = invd * off;
        d_prime[i]      = (b[i] - d_prime[i-1] * off) * invd;
    }

    // last element
    {
        complex_t denom = dp2 - c_prime[N-2] * off;
        complex_t invd  = complex_div(C_ONE, denom);
        d_prime[N-1]    = (b[N-1] - d_prime[N-2] * off) * invd;
    }

    // backward substitution
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}

static void adi_x(const complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N], b[N], x[N];

    complex_t ung   = complex_t(0, dz / (4 * k * dx * dx));
    complex_t neg2u = complex_t(-2.0f, 0.0f) * ung;
    complex_t pos2u = complex_t( 2.0f, 0.0f) * ung;

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            x0[i] = in[i][j];
        }

        complex_t ratio0 = is_near_zero(x0[1])   ? C_ONE : complex_div(x0[0], x0[1]);
        complex_t ratioN = is_near_zero(x0[N-2]) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

        // step B
        complex_t dp_B  = neg2u + C_ONE;
        complex_t dp1_B = dp_B + ung * ratio0;
        complex_t dp2_B = dp_B + ung * ratioN;
        compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        // step A
        complex_t dp_A  = pos2u + C_ONE;
        complex_t dp1_A = dp_A - ung * ratio0;
        complex_t dp2_A = dp_A - ung * ratioN;
        thomas_solver(dp_A, dp1_A, dp2_A, -ung, b, x);

        for (int i = 0; i < N; i++) {
            out[i][j] = x[i];
        }
    }
}

static void adi_y(const complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N], b[N], x[N];

    complex_t ung   = complex_t(0, dz / (4 * k * dy * dy));
    complex_t neg2u = complex_t(-2.0f, 0.0f) * ung;
    complex_t pos2u = complex_t( 2.0f, 0.0f) * ung;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            x0[j] = in[i][j];
        }

        complex_t ratio0 = is_near_zero(x0[1])   ? C_ONE : complex_div(x0[0], x0[1]);
        complex_t ratioN = is_near_zero(x0[N-2]) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

        // step B
        complex_t dp_B  = neg2u + C_ONE;
        complex_t dp1_B = dp_B + ung * ratio0;
        complex_t dp2_B = dp_B + ung * ratioN;
        compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        // step A
        complex_t dp_A  = pos2u + C_ONE;
        complex_t dp1_A = dp_A - ung * ratio0;
        complex_t dp2_A = dp_A - ung * ratioN;
        thomas_solver(dp_A, dp1_A, dp2_A, -ung, b, x);

        for (int j = 0; j < N; j++) {
            out[i][j] = x[j];
        }
    }
}

void diffraction_only(
    const complex_t phi_in[N][N],
          complex_t phi_out[N][N]
) {
#pragma HLS INTERFACE s_axilite port=return bundle=control
#pragma HLS INTERFACE m_axi     port=phi_in  offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE s_axilite port=phi_in  bundle=control
#pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE s_axilite port=phi_out bundle=control

    // Force ONLY these big matrices into URAM; everything else left to inference
    complex_t phi[N][N];
    complex_t tmp[N][N];
#pragma HLS BIND_STORAGE variable=phi type=ram_2p impl=uram
#pragma HLS BIND_STORAGE variable=tmp type=ram_2p impl=uram

    // copy input
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        phi[i][j] = phi_in[i][j];

    // ADI steps (no pipelining)
    adi_x(phi, tmp);
    adi_y(tmp, phi);

    // copy output
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        phi_out[i][j] = phi[i][j];
}
