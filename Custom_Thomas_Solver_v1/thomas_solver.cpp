#include "thomas_solver.h"

void thomas_solver(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t off,
    complex_t b[N],
    complex_t x[N]
) {
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL

#pragma HLS INTERFACE s_axilite port=dp     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=dp1    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=dp2    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=off    bundle=CTRL

#pragma HLS INTERFACE m_axi   port=b   offset=slave   bundle=gmem   depth=N
#pragma HLS INTERFACE mode=s_axilite port=b bundle=control

#pragma HLS INTERFACE m_axi   port=x   offset=slave   bundle=gmem   depth=N
#pragma HLS INTERFACE mode=s_axilite port=x bundle=control

    // buffers en BRAM
    complex_t local_b[N];
#pragma HLS bind_storage variable=local_b type=ram_2p impl=bram
    complex_t c_prime[N];
#pragma HLS bind_storage variable=c_prime type=ram_2p impl=bram
    complex_t d_prime[N];
#pragma HLS bind_storage variable=d_prime type=ram_2p impl=bram

    // 1) carga
    for (int i = 0; i < N; i++) {
        local_b[i] = b[i];
    }

    // 2) forward
    {
        complex_t inv = complex_t(1,0) / dp1;
        c_prime[0] = off * inv;
        d_prime[0] = local_b[0] * inv;
    }
    for (int i = 1; i < N-1; i++) {
        complex_t prev_c = c_prime[i-1];
        complex_t prev_d = d_prime[i-1];
        complex_t denom  = dp - off * prev_c;
        complex_t inv    = complex_t(1,0) / denom;
        c_prime[i] = off * inv;
        d_prime[i] = (local_b[i] - off * prev_d) * inv;
    }
    {
        complex_t denom = dp2 - off * c_prime[N-2];
        complex_t inv   = complex_t(1,0) / denom;
        d_prime[N-1]    = (local_b[N-1] - off * d_prime[N-2]) * inv;
    }

    // 3) backward
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}
