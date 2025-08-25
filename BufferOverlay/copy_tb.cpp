#include "copy_kernel.h"
#include <iostream>

int main() {
  complex_t in[TOTAL];
  complex_t out[TOTAL];

  for (int i = 0; i < TOTAL; ++i) {
    in[i] = complex_t((float)i, (float)(-i));
    out[i] = complex_t(0.0f, 0.0f);
  }

  copy_kernel(in, out);

  for (int i = 0; i < TOTAL; ++i) {
    if (out[i] != in[i]) {
      std::cerr << "Mismatch at index " << i << "\n";
      return 1;
    }
  }
  std::cout << "Copy kernel test passed" << std::endl;
  return 0;
}
