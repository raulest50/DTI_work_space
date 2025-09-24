#include "diff_abs_nsteps.h"
#include <fstream>
#include <iostream>
#include <cmath>

static void read_matrix(const char* fn, complex_t m[N][N]) {
    std::ifstream in(fn);
    if (!in) {
        std::cerr << "Cannot open " << fn << "\n";
        std::exit(1);
    }
    data_t re, im;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            in >> re >> im;
            m[i][j] = complex_t(re, im);
        }
    }
}

int main() {
    complex_t in[N][N];
    complex_t out[N][N];
    complex_t golden[N][N];

    read_matrix("./phi_in.dat", in);
    read_matrix("./golden.dat", golden);

    const data_t alpha = 300.0f;      // 0.3 mm^-1 -> 300 m^-1
    const data_t beta  = 1.0e-11f;    // 10^-11 m/W
    diff_abs_nsteps(in, out, 361, alpha, beta);

    double err = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dre = out[i][j].real() - golden[i][j].real();
            double dim = out[i][j].imag() - golden[i][j].imag();
            err += dre*dre + dim*dim;
        }
    }
    err = std::sqrt(err/(N*N));
    std::cout << "RMS error: " << err << std::endl;
    if (err > 1e-3)
        return 1;
    return 0;
}
