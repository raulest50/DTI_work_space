#include "complex_mul.h"
#include <iostream>

int main() {
    data_t cr = 0.0f;
    data_t ci = 0.0f;
    complex_mul(3.0f, 2.0f, -5.0f, 4.0f, cr, ci);
    const data_t expected_r = -23.0f;
    const data_t expected_i = 2.0f;
    const data_t tol = 1e-5f;
    bool ok = (std::fabs(cr - expected_r) <= tol) && (std::fabs(ci - expected_i) <= tol);
    if (!ok) {
        std::cerr << "Mismatch: got (" << cr << ", " << ci << ")";
        std::cerr << " expected (" << expected_r << ", " << expected_i << ")\n";
        return 1;
    }
    return 0;
}
