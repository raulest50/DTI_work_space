// diffraction_only.cpp
#include "diffraction_only.h"

// Physical constants
static const data_t Lx   = 45e-6f;
static const data_t Ly   = 45e-6f;
static const data_t dx   = Lx / N;
static const data_t dy   = Ly / N;
static const data_t dz   = 1e-4f;
static const data_t k    = 7853981.6339f;
static const data_t eps  = 1e-12f;
static const data_t eps2 = eps * eps;
static const complex_t C_ONE{1.0f, 0.0f};

// “Near zero” check without sqrt or hypot
inline bool is_near_zero(const complex_t &z) {
    data_t re = z.real();
    data_t im = z.imag();
    return (re*re + im*im) < eps2;
}

// Helper: complex divide via conjugate, fully inlined
#pragma HLS INLINE
complex_t complex_div(complex_t num, complex_t den) {
    // compute conjugate of den
    complex_t cd(den.real(), -den.imag());
    // |den|^2
    data_t mag2 = den.real()*den.real() + den.imag()*den.imag();
    // reciprocal
    data_t inv_mag2 = 1.0f / mag2;
    // num * conj(den)
    complex_t p = num * cd;
    // scale result
    return complex_t(p.real() * inv_mag2,
                     p.imag() * inv_mag2);
}

void compute_b_vector(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t off,
    complex_t x0[N],
    complex_t b[N]
) {
    b[0] = dp1 * x0[0] + off * x0[1];
    for (int i = 1; i < N-1; i++) {
        b[i] = off * x0[i-1]
             + dp  * x0[i]
             + off * x0[i+1];
    }
    b[N-1] = off * x0[N-2] + dp2 * x0[N-1];
}

void thomas_solver(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t off,
    complex_t b[N],
    complex_t x[N]
) {
    complex_t c_prime[N];
#pragma HLS bind_storage variable=c_prime type=ram_2p impl=bram
    complex_t d_prime[N];
#pragma HLS bind_storage variable=d_prime type=ram_2p impl=bram

    // first iteration
    complex_t inv1 = complex_div(C_ONE, dp1);
    c_prime[0]    = off * inv1;
    d_prime[0]    = b[0] * inv1;

    // forward sweep
    for (int i = 1; i < N-1; i++) {
        complex_t denom = dp - off * c_prime[i-1];
        complex_t invd  = complex_div(C_ONE, denom);
        c_prime[i]      = off * invd;
        d_prime[i]      = (b[i] - off * d_prime[i-1]) * invd;
    }
    // last element
    {
        complex_t denom = dp2 - off * c_prime[N-2];
        complex_t invd  = complex_div(C_ONE, denom);
        d_prime[N-1]    = (b[N-1] - off * d_prime[N-2]) * invd;
    }

    // backward substitution
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}

void adi_x(complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N], b[N], x[N];
#pragma HLS bind_storage variable=x0 type=ram_2p impl=bram
#pragma HLS bind_storage variable=b   type=ram_2p impl=bram
#pragma HLS bind_storage variable=x   type=ram_2p impl=bram

    // precompute
    complex_t ung   = complex_t(0, dz / (4 * k * dx * dx));
    complex_t neg2u = ung * complex_t(-2.0f, 0.0f);
    complex_t pos2u = ung * complex_t( 2.0f, 0.0f);

    for (int j = 0; j < N; j++) {
        // load column into 1D buffer
        for (int i = 0; i < N; i++) {
            x0[i] = in[i][j];
        }
        // boundary ratios
        complex_t ratio0 = is_near_zero(x0[1])   ? C_ONE : complex_div(x0[0], x0[1]);
        complex_t ratioN = is_near_zero(x0[N-2]) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

        // step B
        complex_t dp_B  = C_ONE + neg2u;
        complex_t dp1_B = dp_B + ung * ratio0;
        complex_t dp2_B = dp_B + ung * ratioN;
        compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        // step A
        complex_t dp_A  = C_ONE + pos2u;
        complex_t dp1_A = dp_A - ung * ratio0;
        complex_t dp2_A = dp_A - ung * ratioN;
        thomas_solver(dp_A, dp1_A, dp2_A, -ung, b, x);

        // store result
        for (int i = 0; i < N; i++) {
            out[i][j] = x[i];
        }
    }
}

void adi_y(complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N], b[N], x[N];
#pragma HLS bind_storage variable=x0 type=ram_2p impl=bram
#pragma HLS bind_storage variable=b   type=ram_2p impl=bram
#pragma HLS bind_storage variable=x   type=ram_2p impl=bram

    // precompute
    complex_t ung   = complex_t(0, dz / (4 * k * dy * dy));
    complex_t neg2u = ung * complex_t(-2.0f, 0.0f);
    complex_t pos2u = ung * complex_t( 2.0f, 0.0f);

    for (int i = 0; i < N; i++) {
        // load row into 1D buffer
        for (int j = 0; j < N; j++) {
            x0[j] = in[i][j];
        }
        // boundary ratios
        complex_t ratio0 = is_near_zero(x0[1])   ? C_ONE : complex_div(x0[0], x0[1]);
        complex_t ratioN = is_near_zero(x0[N-2]) ? C_ONE : complex_div(x0[N-1], x0[N-2]);

        // step B
        complex_t dp_B  = C_ONE + neg2u;
        complex_t dp1_B = dp_B + ung * ratio0;
        complex_t dp2_B = dp_B + ung * ratioN;
        compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        // step A
        complex_t dp_A  = C_ONE + pos2u;
        complex_t dp1_A = dp_A - ung * ratio0;
        complex_t dp2_A = dp_A - ung * ratioN;
        thomas_solver(dp_A, dp1_A, dp2_A, -ung, b, x);

        // store result
        for (int j = 0; j < N; j++) {
            out[i][j] = x[j];
        }
    }
}

void diffraction_only(
    complex_t phi_in[N][N],
    complex_t phi_out[N][N]
) {
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE m_axi     port=phi_in  offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE s_axilite port=phi_in  bundle=control
#pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE s_axilite port=phi_out bundle=control

    complex_t phi[N][N];
    complex_t tmp[N][N];
#pragma HLS bind_storage variable=phi type=ram_2p impl=uram
#pragma HLS bind_storage variable=tmp type=ram_2p impl=uram

    // copy input
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        phi[i][j] = phi_in[i][j];

    // ADI steps
    adi_x(phi, tmp);
    adi_y(tmp, phi);

    // copy output
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        phi_out[i][j] = phi[i][j];
}
