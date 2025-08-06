#include "thomas_solver.h"
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>

static void read_vec(const char* fn, std::vector<complex_t>& v) {
    std::ifstream in(fn);
    if (!in) {
        std::cerr << "Cannot open " << fn << "\n";
        std::exit(1);
    }
    data_t re, im;
    while (in >> re >> im)
        v.emplace_back(re, im);
}

static void write_vec(const char* fn, const std::vector<complex_t>& v) {
    std::ofstream out(fn);
    for (auto& z : v)
        out << z.real() << ' ' << z.imag() << '\n';
}

// CPU reference using std::complex<float>
static std::vector<std::complex<float>> ref_solver(
        std::complex<float> dp,
        std::complex<float> dp1,
        std::complex<float> dp2,
        std::complex<float> off,
        const std::vector<std::complex<float>>& b) {
    int n = b.size();
    std::vector<std::complex<float>> cprime(n);
    std::vector<std::complex<float>> dprime(n);
    std::vector<std::complex<float>> x(n);

    cprime[0] = off / dp1;
    dprime[0] = b[0] / dp1;
    for (int i=1;i<n-1;i++) {
        std::complex<float> denom = dp - off*cprime[i-1];
        cprime[i] = off/denom;
        dprime[i] = (b[i]-off*dprime[i-1])/denom;
    }
    dprime[n-1] = (b[n-1] - off*dprime[n-2]) /
                  (dp2 - off*cprime[n-2]);

    x[n-1]=dprime[n-1];
    for(int i=n-2;i>=0;i--) {
        x[i]=dprime[i]-cprime[i]*x[i+1];
    }
    return x;
}

int main() {
    // Constant tridiagonal coefficients
    complex_t dp  = complex_t(5.0f, 0.5f);
    complex_t dp1 = complex_t(5.5f, 0.5f);
    complex_t dp2 = complex_t(4.5f, 0.5f);
    complex_t off = complex_t(1.0f, -0.2f);
    std::vector<complex_t> B;
    read_vec("./b.dat", B);
    assert(B.size()==N);

    // Run HLS function
    std::vector<complex_t> X(N);
    thomas_solver(dp, dp1, dp2, off, B.data(), X.data());

    // Compute reference on CPU for comparison
    std::vector<std::complex<float>> b_f(N);
    for(int i=0;i<N;i++) {
        b_f[i] = std::complex<float>(float(B[i].real()), float(B[i].imag()));
    }
    auto ref = ref_solver(
            std::complex<float>(float(dp.real()), float(dp.imag())),
            std::complex<float>(float(dp1.real()), float(dp1.imag())),
            std::complex<float>(float(dp2.real()), float(dp2.imag())),
            std::complex<float>(float(off.real()), float(off.imag())),
            b_f);

    double err=0.0;
    for(int i=0;i<N;i++) {
        double dre = float(X[i].real()) - ref[i].real();
        double dim = float(X[i].imag()) - ref[i].imag();
        err += dre*dre + dim*dim;
    }
    err = std::sqrt(err/N);
    std::cout << "RMS error: " << err << std::endl;

    std::vector<complex_t> ref_out(N);
    for(int i=0;i<N;i++)
        ref_out[i] = complex_t(data_t(ref[i].real()), data_t(ref[i].imag()));
    write_vec("./golden.dat", ref_out);
    return 0;
}
