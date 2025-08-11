#include "complex_div.h"
#include <iostream>
#include <cmath>

int main() {
    complex_t a(4.0f, 2.0f);
    complex_t b(1.0f, -3.0f);
    complex_t c;

    complex_div(a, b, c);

    complex_t expected = a / b;
    float err_r = expected.real() - c.real();
    float err_i = expected.imag() - c.imag();
    float rmse = std::sqrt((err_r * err_r + err_i * err_i) / 2.0f);

    std::cout << "RMSE: " << rmse << "\n";
    const float tol = 1e-5f;
    if (rmse > tol) {
        std::cerr << "RMS error too high\n";
        return 1;
    }
    return 0;
}
