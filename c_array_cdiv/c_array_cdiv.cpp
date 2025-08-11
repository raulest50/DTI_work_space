#include "c_array_cdiv.h"

static complex_t complex_div(const complex_t &a, const complex_t &b) {
    data_t br = b.real();
    data_t bi = b.imag();
    data_t denom = br*br + bi*bi;
    data_t inv   = (data_t)1.0f / denom;
    complex_t conj_b(br, -bi);
    complex_t num = a * conj_b;
    return complex_t(num.real() * inv, num.imag() * inv);
}

static void array_to_stream(const complex_t in[N], hls::stream<complex_t> &s) {
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        s.write(in[i]);
    }
}

static void stream_div(hls::stream<complex_t> &in_s,
                       hls::stream<complex_t> &out_s,
                       const complex_t &div) {
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        complex_t v = in_s.read();
        out_s.write(complex_div(v, div));
    }
}

static void stream_to_array(hls::stream<complex_t> &s, complex_t out[N]) {
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        out[i] = s.read();
    }
}

void c_array_cdiv(const complex_t in[N],
                  complex_t out[N],
                  data_t div_real,
                  data_t div_imag) {
#pragma HLS INTERFACE s_axilite port=return bundle=control
#pragma HLS INTERFACE m_axi     port=in  offset=slave bundle=gmem depth=64
#pragma HLS INTERFACE s_axilite port=in  bundle=control
#pragma HLS INTERFACE m_axi     port=out offset=slave bundle=gmem depth=64
#pragma HLS INTERFACE s_axilite port=out bundle=control
#pragma HLS INTERFACE s_axilite port=div_real bundle=control
#pragma HLS INTERFACE s_axilite port=div_imag bundle=control

    hls::stream<complex_t> s_in;
    hls::stream<complex_t> s_out;
#pragma HLS STREAM variable=s_in depth=64
#pragma HLS STREAM variable=s_out depth=64

    complex_t div(div_real, div_imag);

    array_to_stream(in, s_in);
    stream_div(s_in, s_out, div);
    stream_to_array(s_out, out);
}
