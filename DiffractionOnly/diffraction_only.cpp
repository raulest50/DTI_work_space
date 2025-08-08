#include "diffraction_only.h"

// Physical constants
static const data_t Lx = 45e-6f;
static const data_t Ly = 45e-6f;
static const data_t dx = Lx / N;
static const data_t dy = Ly / N;
static const data_t dz = 1e-4f;
static const data_t k  = 7853981.6339f;
static const data_t eps = 1e-12f;

inline data_t abs_complex(complex_t z) {
    return hls::hypot(z.real(), z.imag());
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
    for (int i = 1; i < N - 1; i++) {
        b[i] = off * x0[i-1] + dp * x0[i] + off * x0[i+1];
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

    complex_t inv = complex_t(1,0) / dp1;
    c_prime[0] = off * inv;
    d_prime[0] = b[0] * inv;

    for (int i = 1; i < N-1; i++) {
        complex_t denom = dp - off * c_prime[i-1];
        complex_t inv_d = complex_t(1,0) / denom;
        c_prime[i] = off * inv_d;
        d_prime[i] = (b[i] - off * d_prime[i-1]) * inv_d;
    }

    complex_t denom = dp2 - off * c_prime[N-2];
    complex_t inv_d = complex_t(1,0) / denom;
    d_prime[N-1] = (b[N-1] - off * d_prime[N-2]) * inv_d;

    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}

void adi_x(complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N];
    complex_t b[N];
    complex_t x[N];
#pragma HLS bind_storage variable=x0 type=ram_2p impl=bram
#pragma HLS bind_storage variable=b  type=ram_2p impl=bram
#pragma HLS bind_storage variable=x  type=ram_2p impl=bram
    complex_t ung(0, dz / (4 * k * dx * dx));
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            x0[i] = in[i][j];
        }
        complex_t ratio_x0 = abs_complex(x0[1]) < eps ? complex_t(1,0) : x0[0] / x0[1];
        complex_t ratio_xn = abs_complex(x0[N-2]) < eps ? complex_t(1,0) : x0[N-1] / x0[N-2];

        complex_t dp1_B = -2.0f*ung + complex_t(1,0) + ung * ratio_x0;
        complex_t dp2_B = -2.0f*ung + complex_t(1,0) + ung * ratio_xn;
        complex_t dp_B  = -2.0f*ung + complex_t(1,0);
        complex_t do_B  = ung;

        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, x0, b);

        complex_t dp1_A = 2.0f*ung + complex_t(1,0) - ung * ratio_x0;
        complex_t dp2_A = 2.0f*ung + complex_t(1,0) - ung * ratio_xn;
        complex_t dp_A  = 2.0f*ung + complex_t(1,0);
        complex_t do_A  = -ung;

        thomas_solver(dp_A, dp1_A, dp2_A, do_A, b, x);

        for (int i = 0; i < N; i++) {
            out[i][j] = x[i];
        }
    }
}

void adi_y(complex_t in[N][N], complex_t out[N][N]) {
    complex_t x0[N];
    complex_t b[N];
    complex_t x[N];
#pragma HLS bind_storage variable=x0 type=ram_2p impl=bram
#pragma HLS bind_storage variable=b  type=ram_2p impl=bram
#pragma HLS bind_storage variable=x  type=ram_2p impl=bram
    complex_t ung(0, dz / (4 * k * dy * dy));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            x0[j] = in[i][j];
        }
        complex_t ratio_y0 = abs_complex(x0[1]) < eps ? complex_t(1,0) : x0[0] / x0[1];
        complex_t ratio_yn = abs_complex(x0[N-2]) < eps ? complex_t(1,0) : x0[N-1] / x0[N-2];

        complex_t dp1_B = -2.0f*ung + complex_t(1,0) + ung * ratio_y0;
        complex_t dp2_B = -2.0f*ung + complex_t(1,0) + ung * ratio_yn;
        complex_t dp_B  = -2.0f*ung + complex_t(1,0);
        complex_t do_B  = ung;

        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, x0, b);

        complex_t dp1_A = 2.0f*ung + complex_t(1,0) - ung * ratio_y0;
        complex_t dp2_A = 2.0f*ung + complex_t(1,0) - ung * ratio_yn;
        complex_t dp_A  = 2.0f*ung + complex_t(1,0);
        complex_t do_A  = -ung;

        thomas_solver(dp_A, dp1_A, dp2_A, do_A, b, x);

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
#pragma HLS INTERFACE m_axi   port=phi_in  offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE mode=s_axilite port=phi_in  bundle=control
#pragma HLS INTERFACE m_axi   port=phi_out offset=slave bundle=gmem depth=4096
#pragma HLS INTERFACE mode=s_axilite port=phi_out bundle=control

    complex_t phi[N][N];
    complex_t tmp[N][N];
#pragma HLS bind_storage variable=phi type=ram_2p impl=uram
#pragma HLS bind_storage variable=tmp type=ram_2p impl=uram

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            phi[i][j] = phi_in[i][j];
        }
    }

    adi_x(phi, tmp);
    adi_y(tmp, phi);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            phi_out[i][j] = phi[i][j];
        }
    }
}

