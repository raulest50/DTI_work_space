#include "copy_kernel.h"

extern "C" void copy_kernel(const complex_t* in, complex_t* out) {
#pragma HLS INTERFACE s_axilite port=return bundle=control
#pragma HLS INTERFACE m_axi port=in  bundle=gmem depth=4096
#pragma HLS INTERFACE m_axi port=out bundle=gmem depth=4096
#pragma HLS INTERFACE s_axilite port=in  bundle=control
#pragma HLS INTERFACE s_axilite port=out bundle=control

  for (int idx = 0; idx < TOTAL; ++idx) {
#pragma HLS PIPELINE II=1
    out[idx] = in[idx];
  }
}
