#include "diffraction_only.h"

int main() {
    complex_t in[N][N];
    complex_t out[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            in[i][j] = complex_t(0, 0);
        }
    }
    diffraction_only(in, out);
    return 0;
}
