// diff_losses_tb.cpp - Diagnostic testbench comparing float vs fixed-point
#include "diff_losses.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <complex>

// ============================================================================
// FLOAT REFERENCE IMPLEMENTATION (for comparison)
// ============================================================================
using float_t = float;
using complex_f = std::complex<float>;

static const float_t F_UNG_COEFF = 6.4368f;
static const float_t F_DZ_OVER_4 = 2.5e-5f;
static const float_t F_EPS2 = 1e-12f;
static const complex_f F_ONE = complex_f(1.0f, 0.0f);
static const complex_f F_ZERO = complex_f(0.0f, 0.0f);

// Float wrappers
static complex_f f_add(const complex_f &a, const complex_f &b) { return a + b; }
static complex_f f_sub(const complex_f &a, const complex_f &b) { return a - b; }
static complex_f f_mul(const complex_f &a, const complex_f &b) { return a * b; }
static complex_f f_scale(const complex_f &a, float_t s) { return complex_f(a.real()*s, a.imag()*s); }
static float_t f_abs2(const complex_f &a) { return a.real()*a.real() + a.imag()*a.imag(); }
static complex_f f_conj(const complex_f &a) { return std::conj(a); }
static bool f_near_zero(const complex_f &z, float_t eps2) { return f_abs2(z) < eps2; }

static complex_f f_div_safe(const complex_f &num, const complex_f &den, float_t eps2) {
    float_t mag2 = f_abs2(den);
    if (mag2 < eps2) mag2 = eps2;
    float_t inv = 1.0f / mag2;
    complex_f p = f_mul(f_conj(den), num);
    return complex_f(p.real()*inv, p.imag()*inv);
}

static void f_half_linear_absorption(complex_f phi[N][N], float_t alpha) {
    float_t c = std::exp(-alpha * F_DZ_OVER_4);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            phi[i][j] = f_scale(phi[i][j], c);
}

static void f_half_2photon_absorption(complex_f phi[N][N], float_t beta) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            float_t atten = std::exp(-beta * F_DZ_OVER_4 * f_abs2(phi[i][j]));
            phi[i][j] = f_scale(phi[i][j], atten);
        }
}

static void f_apply_N_midpoint(complex_f phi[N][N], float_t alpha, float_t beta) {
    f_half_linear_absorption(phi, alpha);
    f_half_2photon_absorption(phi, beta);
    f_half_2photon_absorption(phi, beta);
    f_half_linear_absorption(phi, alpha);
}

static void f_compute_b_vector(complex_f dp, complex_f dp1, complex_f dp2, complex_f off,
                               const complex_f x0[N], complex_f b[N]) {
    b[0] = f_add(f_mul(x0[0], dp1), f_mul(x0[1], off));
    for (int i = 1; i < N-1; ++i) {
        complex_f t = F_ZERO;
        t = f_add(t, f_mul(x0[i-1], off));
        t = f_add(t, f_mul(x0[i], dp));
        t = f_add(t, f_mul(x0[i+1], off));
        b[i] = t;
    }
    b[N-1] = f_add(f_mul(x0[N-2], off), f_mul(x0[N-1], dp2));
}

static void f_thomas_solver(complex_f dp, complex_f dp1, complex_f dp2, complex_f off,
                            complex_f b[N], complex_f x[N]) {
    complex_f c_prime[N], d_prime[N];

    complex_f inv1 = f_div_safe(F_ONE, dp1, F_EPS2);
    c_prime[0] = f_mul(inv1, off);
    d_prime[0] = f_mul(b[0], inv1);

    for (int i = 1; i < N-1; ++i) {
        complex_f denom = f_sub(dp, f_mul(c_prime[i-1], off));
        complex_f invd = f_div_safe(F_ONE, denom, F_EPS2);
        c_prime[i] = f_mul(invd, off);
        d_prime[i] = f_mul(f_sub(b[i], f_mul(d_prime[i-1], off)), invd);
    }

    complex_f denom = f_sub(dp2, f_mul(c_prime[N-2], off));
    complex_f invd = f_div_safe(F_ONE, denom, F_EPS2);
    d_prime[N-1] = f_mul(f_sub(b[N-1], f_mul(d_prime[N-2], off)), invd);

    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; --i)
        x[i] = f_sub(d_prime[i], f_mul(c_prime[i], x[i+1]));
}

static void f_adi_x(const complex_f in[N][N], complex_f out[N][N]) {
    complex_f x0[N], b[N], x[N];
    complex_f ung = complex_f(0.0f, F_UNG_COEFF);
    complex_f neg2u = f_mul(complex_f(-2.0f, 0.0f), ung);
    complex_f pos2u = f_mul(complex_f(2.0f, 0.0f), ung);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) x0[i] = in[i][j];

        complex_f ratio0 = f_near_zero(x0[1], F_EPS2) ? F_ONE : f_div_safe(x0[0], x0[1], F_EPS2);
        complex_f ratioN = f_near_zero(x0[N-2], F_EPS2) ? F_ONE : f_div_safe(x0[N-1], x0[N-2], F_EPS2);

        complex_f dp_B = f_add(neg2u, F_ONE);
        complex_f dp1_B = f_add(dp_B, f_mul(ung, ratio0));
        complex_f dp2_B = f_add(dp_B, f_mul(ung, ratioN));
        f_compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        complex_f dp_A = f_add(pos2u, F_ONE);
        complex_f dp1_A = f_sub(dp_A, f_mul(ung, ratio0));
        complex_f dp2_A = f_sub(dp_A, f_mul(ung, ratioN));
        f_thomas_solver(dp_A, dp1_A, dp2_A, f_scale(ung, -1.0f), b, x);

        for (int i = 0; i < N; ++i) out[i][j] = x[i];
    }
}

static void f_adi_y(const complex_f in[N][N], complex_f out[N][N]) {
    complex_f x0[N], b[N], x[N];
    complex_f ung = complex_f(0.0f, F_UNG_COEFF);
    complex_f neg2u = f_mul(complex_f(-2.0f, 0.0f), ung);
    complex_f pos2u = f_mul(complex_f(2.0f, 0.0f), ung);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) x0[j] = in[i][j];

        complex_f ratio0 = f_near_zero(x0[1], F_EPS2) ? F_ONE : f_div_safe(x0[0], x0[1], F_EPS2);
        complex_f ratioN = f_near_zero(x0[N-2], F_EPS2) ? F_ONE : f_div_safe(x0[N-1], x0[N-2], F_EPS2);

        complex_f dp_B = f_add(neg2u, F_ONE);
        complex_f dp1_B = f_add(dp_B, f_mul(ung, ratio0));
        complex_f dp2_B = f_add(dp_B, f_mul(ung, ratioN));
        f_compute_b_vector(dp_B, dp1_B, dp2_B, ung, x0, b);

        complex_f dp_A = f_add(pos2u, F_ONE);
        complex_f dp1_A = f_sub(dp_A, f_mul(ung, ratio0));
        complex_f dp2_A = f_sub(dp_A, f_mul(ung, ratioN));
        f_thomas_solver(dp_A, dp1_A, dp2_A, f_scale(ung, -1.0f), b, x);

        for (int j = 0; j < N; ++j) out[i][j] = x[j];
    }
}

// ============================================================================
// FIXED-POINT WRAPPERS (duplicated from diff_losses.cpp since they're static)
// ============================================================================

static complex_t c_add(const complex_t &a, const complex_t &b) {
    return complex_t(a.real()+b.real(), a.imag()+b.imag());
}
static complex_t c_sub(const complex_t &a, const complex_t &b) {
    return complex_t(a.real()-b.real(), a.imag()-b.imag());
}
static complex_t c_scale(const complex_t &a, data_t s) {
    return complex_t(a.real()*s, a.imag()*s);
}
static complex_t c_mul(const complex_t &a, const complex_t &b) {
    data_t ar = a.real(), ai = a.imag();
    data_t br = b.real(), bi = b.imag();
    return complex_t(ar*br - ai*bi, ar*bi + ai*br);
}
static data_t c_abs2(const complex_t &a) {
    data_t ar = a.real(), ai = a.imag();
    return ar*ar + ai*ai;
}
static complex_t c_conj(const complex_t &a) {
    return complex_t(a.real(), -a.imag());
}
static complex_t c_div_safe(const complex_t &num, const complex_t &den, data_t eps2) {
    data_t mag2 = c_abs2(den);
    if (mag2 < eps2) mag2 = eps2;
    data_t inv = data_t(1) / mag2;
    complex_t p = c_mul(c_conj(den), num);
    return complex_t(p.real()*inv, p.imag()*inv);
}
static bool c_near_zero(const complex_t &z, data_t eps2) {
    return c_abs2(z) < eps2;
}

// ============================================================================
// COMPARISON UTILITIES
// ============================================================================

static double to_double(data_t x) { return (double)x; }

static void read_matrix_float(const char* fn, complex_f m[N][N]) {
    std::ifstream in(fn);
    if (!in) { std::cerr << "Cannot open " << fn << "\n"; std::exit(1); }
    float re, im;
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++) {
            in >> re >> im;
            m[i][j] = complex_f(re, im);
        }
}

static void read_matrix_fixed(const char* fn, complex_t m[N][N]) {
    std::ifstream in(fn);
    if (!in) { std::cerr << "Cannot open " << fn << "\n"; std::exit(1); }
    double re, im;
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++) {
            in >> re >> im;
            m[i][j] = complex_t(data_t(re), data_t(im));
        }
}

static void copy_float_to_fixed(const complex_f src[N][N], complex_t dst[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            dst[i][j] = complex_t(data_t(src[i][j].real()), data_t(src[i][j].imag()));
}

static void copy_fixed_to_float(const complex_t src[N][N], complex_f dst[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            dst[i][j] = complex_f((float)src[i][j].real(), (float)src[i][j].imag());
}

// Compute RMS difference between float and fixed arrays
static double compute_rms_diff(const complex_f flt[N][N], const complex_t fix[N][N]) {
    double err = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dre = flt[i][j].real() - to_double(fix[i][j].real());
            double dim = flt[i][j].imag() - to_double(fix[i][j].imag());
            err += dre*dre + dim*dim;
        }
    }
    return std::sqrt(err / (N*N));
}

// Log sample values from center of grid
static void log_samples(std::ofstream& log, const char* label,
                        const complex_f flt[N][N], const complex_t fix[N][N]) {
    log << "\n" << label << ":\n";
    log << std::scientific << std::setprecision(6);

    // Log center pixel and a few neighbors
    int ci = N/2, cj = N/2;
    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            int i = ci + di, j = cj + dj;
            log << "  [" << i << "," << j << "] float=("
                << flt[i][j].real() << ", " << flt[i][j].imag() << ") "
                << "fixed=(" << to_double(fix[i][j].real()) << ", "
                << to_double(fix[i][j].imag()) << ")\n";
        }
    }

    double rms = compute_rms_diff(flt, fix);
    log << "  RMS diff: " << rms << "\n";
}

// ============================================================================
// MAIN DIAGNOSTIC TEST
// ============================================================================

int main() {
    std::ofstream log("fp_debug.txt");
    if (!log) { std::cerr << "Cannot open fp_debug.txt\n"; return 1; }

    log << "============================================\n";
    log << "Fixed-Point vs Float Diagnostic Comparison\n";
    log << "============================================\n";
    log << "data_t: ap_fixed<" << DATA_TOTAL_BITS << ", " << DATA_INT_BITS << ">\n";
    log << "Fractional bits: " << (DATA_TOTAL_BITS - DATA_INT_BITS) << "\n";
    log << "Resolution: " << std::pow(2.0, -(DATA_TOTAL_BITS - DATA_INT_BITS)) << "\n";
    log << "============================================\n";

    // Load input data
    complex_f phi_f[N][N], tmp_f[N][N];
    complex_t phi_x[N][N], tmp_x[N][N];

    read_matrix_float("./phi_in.dat", phi_f);
    read_matrix_fixed("./phi_in.dat", phi_x);

    const int TEST_STEPS = 10;
    const float_t alpha = 0.0f;
    const float_t beta = 0.0f;
    const data_t alpha_x = data_t(0);
    const data_t beta_x = data_t(0);

    log << "\nRunning " << TEST_STEPS << " steps with alpha=0, beta=0\n";
    log << "Comparing float vs fixed-point at each stage\n";

    // Log initial state
    log_samples(log, "INITIAL INPUT", phi_f, phi_x);

    for (int s = 0; s < TEST_STEPS; s++) {
        std::cout << "Step " << s << "/" << TEST_STEPS << std::endl;

        bool log_this_step = (s == 0 || s == TEST_STEPS - 1);

        // === ADI X ===
        f_adi_x(phi_f, tmp_f);

        // For fixed: we need to replicate adi_x logic or call internal function
        // Since adi_x is static in diff_losses.cpp, we'll replicate here
        {
            complex_t x0[N], b[N], x[N];
            complex_t ung = c_make(data_t(0), data_t(6.4368));
            complex_t neg2u = c_mul(c_make(data_t(-2), data_t(0)), ung);
            complex_t pos2u = c_mul(c_make(data_t(2), data_t(0)), ung);
            data_t eps2_local = data_t(1e-5);

            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < N; ++i) x0[i] = phi_x[i][j];

                complex_t ratio0 = c_near_zero(x0[1], eps2_local) ? c_make(data_t(1), data_t(0)) :
                    c_div_safe(x0[0], x0[1], eps2_local);
                complex_t ratioN = c_near_zero(x0[N-2], eps2_local) ? c_make(data_t(1), data_t(0)) :
                    c_div_safe(x0[N-1], x0[N-2], eps2_local);

                complex_t dp_B = c_add(neg2u, c_make(data_t(1), data_t(0)));
                complex_t dp1_B = c_add(dp_B, c_mul(ung, ratio0));
                complex_t dp2_B = c_add(dp_B, c_mul(ung, ratioN));

                // compute b vector
                b[0] = c_add(c_mul(x0[0], dp1_B), c_mul(x0[1], ung));
                for (int i = 1; i < N-1; ++i) {
                    complex_t t = c_make(data_t(0), data_t(0));
                    t = c_add(t, c_mul(x0[i-1], ung));
                    t = c_add(t, c_mul(x0[i], dp_B));
                    t = c_add(t, c_mul(x0[i+1], ung));
                    b[i] = t;
                }
                b[N-1] = c_add(c_mul(x0[N-2], ung), c_mul(x0[N-1], dp2_B));

                complex_t dp_A = c_add(pos2u, c_make(data_t(1), data_t(0)));
                complex_t dp1_A = c_sub(dp_A, c_mul(ung, ratio0));
                complex_t dp2_A = c_sub(dp_A, c_mul(ung, ratioN));
                complex_t off_A = c_scale(ung, data_t(-1));

                // thomas solver
                complex_t c_prime[N], d_prime[N];
                complex_t inv1 = c_div_safe(c_make(data_t(1), data_t(0)), dp1_A, eps2_local);
                c_prime[0] = c_mul(inv1, off_A);
                d_prime[0] = c_mul(b[0], inv1);

                for (int i = 1; i < N-1; ++i) {
                    complex_t denom = c_sub(dp_A, c_mul(c_prime[i-1], off_A));
                    complex_t invd = c_div_safe(c_make(data_t(1), data_t(0)), denom, eps2_local);
                    c_prime[i] = c_mul(invd, off_A);
                    d_prime[i] = c_mul(c_sub(b[i], c_mul(d_prime[i-1], off_A)), invd);
                }

                complex_t denom = c_sub(dp2_A, c_mul(c_prime[N-2], off_A));
                complex_t invd = c_div_safe(c_make(data_t(1), data_t(0)), denom, eps2_local);
                d_prime[N-1] = c_mul(c_sub(b[N-1], c_mul(d_prime[N-2], off_A)), invd);

                x[N-1] = d_prime[N-1];
                for (int i = N-2; i >= 0; --i)
                    x[i] = c_sub(d_prime[i], c_mul(c_prime[i], x[i+1]));

                for (int i = 0; i < N; ++i) tmp_x[i][j] = x[i];
            }
        }

        if (log_this_step) {
            char label[64];
            sprintf(label, "STEP %d - After ADI_X", s);
            log_samples(log, label, tmp_f, tmp_x);
        }

        // === Apply N midpoint (absorption) ===
        f_apply_N_midpoint(tmp_f, alpha, beta);
        // For fixed - with alpha=beta=0, this is a no-op (exp(0)=1)
        // But let's still apply it to match
        {
            data_t c = data_t(expf(-float(alpha_x) * 2.5e-5f));
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    tmp_x[i][j] = c_scale(tmp_x[i][j], c);
            // TPA (no-op when beta=0)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    data_t atten = data_t(expf(-float(beta_x) * 2.5e-5f * float(c_abs2(tmp_x[i][j]))));
                    tmp_x[i][j] = c_scale(tmp_x[i][j], atten);
                }
            // Repeat TPA
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    data_t atten = data_t(expf(-float(beta_x) * 2.5e-5f * float(c_abs2(tmp_x[i][j]))));
                    tmp_x[i][j] = c_scale(tmp_x[i][j], atten);
                }
            // Second linear absorption
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    tmp_x[i][j] = c_scale(tmp_x[i][j], c);
        }

        if (log_this_step) {
            char label[64];
            sprintf(label, "STEP %d - After ABSORPTION", s);
            log_samples(log, label, tmp_f, tmp_x);
        }

        // === ADI Y ===
        f_adi_y(tmp_f, phi_f);

        // Fixed ADI Y
        {
            complex_t x0[N], b[N], x[N];
            complex_t ung = c_make(data_t(0), data_t(6.4368));
            complex_t neg2u = c_mul(c_make(data_t(-2), data_t(0)), ung);
            complex_t pos2u = c_mul(c_make(data_t(2), data_t(0)), ung);
            data_t eps2_local = data_t(1e-5);

            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) x0[j] = tmp_x[i][j];

                complex_t ratio0 = c_near_zero(x0[1], eps2_local) ? c_make(data_t(1), data_t(0)) :
                    c_div_safe(x0[0], x0[1], eps2_local);
                complex_t ratioN = c_near_zero(x0[N-2], eps2_local) ? c_make(data_t(1), data_t(0)) :
                    c_div_safe(x0[N-1], x0[N-2], eps2_local);

                complex_t dp_B = c_add(neg2u, c_make(data_t(1), data_t(0)));
                complex_t dp1_B = c_add(dp_B, c_mul(ung, ratio0));
                complex_t dp2_B = c_add(dp_B, c_mul(ung, ratioN));

                b[0] = c_add(c_mul(x0[0], dp1_B), c_mul(x0[1], ung));
                for (int jj = 1; jj < N-1; ++jj) {
                    complex_t t = c_make(data_t(0), data_t(0));
                    t = c_add(t, c_mul(x0[jj-1], ung));
                    t = c_add(t, c_mul(x0[jj], dp_B));
                    t = c_add(t, c_mul(x0[jj+1], ung));
                    b[jj] = t;
                }
                b[N-1] = c_add(c_mul(x0[N-2], ung), c_mul(x0[N-1], dp2_B));

                complex_t dp_A = c_add(pos2u, c_make(data_t(1), data_t(0)));
                complex_t dp1_A = c_sub(dp_A, c_mul(ung, ratio0));
                complex_t dp2_A = c_sub(dp_A, c_mul(ung, ratioN));
                complex_t off_A = c_scale(ung, data_t(-1));

                complex_t c_prime[N], d_prime[N];
                complex_t inv1 = c_div_safe(c_make(data_t(1), data_t(0)), dp1_A, eps2_local);
                c_prime[0] = c_mul(inv1, off_A);
                d_prime[0] = c_mul(b[0], inv1);

                for (int jj = 1; jj < N-1; ++jj) {
                    complex_t denom = c_sub(dp_A, c_mul(c_prime[jj-1], off_A));
                    complex_t invd = c_div_safe(c_make(data_t(1), data_t(0)), denom, eps2_local);
                    c_prime[jj] = c_mul(invd, off_A);
                    d_prime[jj] = c_mul(c_sub(b[jj], c_mul(d_prime[jj-1], off_A)), invd);
                }

                complex_t denom = c_sub(dp2_A, c_mul(c_prime[N-2], off_A));
                complex_t invd = c_div_safe(c_make(data_t(1), data_t(0)), denom, eps2_local);
                d_prime[N-1] = c_mul(c_sub(b[N-1], c_mul(d_prime[N-2], off_A)), invd);

                x[N-1] = d_prime[N-1];
                for (int jj = N-2; jj >= 0; --jj)
                    x[jj] = c_sub(d_prime[jj], c_mul(c_prime[jj], x[jj+1]));

                for (int j = 0; j < N; ++j) phi_x[i][j] = x[j];
            }
        }

        if (log_this_step) {
            char label[64];
            sprintf(label, "STEP %d - After ADI_Y (end of step)", s);
            log_samples(log, label, phi_f, phi_x);
        }
    }

    // Final comparison
    double final_rms = compute_rms_diff(phi_f, phi_x);

    log << "\n============================================\n";
    log << "FINAL RESULTS after " << TEST_STEPS << " steps:\n";
    log << "RMS difference (float vs fixed): " << std::scientific << final_rms << "\n";
    log << "============================================\n";

    log.close();

    std::cout << "\nDiagnostic complete. Results in fp_debug.txt\n";
    std::cout << "Final RMS diff: " << std::scientific << final_rms << std::endl;

    if (final_rms > 1.0) {
        std::cout << "FAIL: Large divergence detected\n";
        return 1;
    }

    std::cout << "PASS: Float and fixed-point are reasonably close\n";
    return 0;
}
