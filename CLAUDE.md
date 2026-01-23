# CLAUDE.md - Vitis HLS Project

---
---
# PART 1: GENERIC VITIS HLS DIRECTIVES (Reusable)
---
---

## Tool Context

This is a **Xilinx/AMD Vitis HLS 2025.2** project. Code written here is **C++ for hardware synthesis**, NOT standard C++. The HLS compiler transforms C++ into RTL (Verilog/VHDL) for FPGA implementation.

---

## HLS C++ Constraints (CRITICAL)

### Forbidden Constructs - NEVER USE
| Construct | Reason |
|-----------|--------|
| **Recursion** | Hardware cannot implement recursive call stacks |
| **Dynamic memory** (`new`, `delete`, `malloc`, `free`) | No heap in hardware |
| **System calls** (`printf`, `scanf`, file I/O) | No OS in hardware (OK in testbench only) |
| **Virtual functions / polymorphism** | Cannot resolve at synthesis time |
| **Exceptions** (`try`/`catch`/`throw`) | No exception handling in hardware |
| **STL containers** (`std::vector`, `std::map`, `std::string`) | Dynamic allocation |
| **Function pointers** | Limited/no support |
| **Variable-length arrays** | Array sizes must be compile-time constants |

### Required Patterns
- Use **fixed-size arrays** with compile-time constants: `data_t arr[N][M]`
- Use **`constexpr`** or **`#define`** for array dimensions
- Use **static functions** with `#pragma HLS INLINE off` to control resource sharing
- All loop bounds should be determinable at compile time (or use `LOOP_TRIPCOUNT`)

---

## Complex Numbers in Vitis HLS

### Summary Table

| Type | Simulation | Synthesis | Recommendation |
|------|------------|-----------|----------------|
| `std::complex<ap_fixed<W,I>>` | Yes | Yes | **Best for fixed-point** |
| `std::complex<float>` | Yes | **Has issues** | Use wrappers or avoid |
| `std::complex<double>` | Yes | **Has issues** | Use wrappers or avoid |
| `complex_wrapper<float>` | Yes | Yes | **Best for float** (Vitis DSP Lib) |
| `complex_wrapper<double>` | Yes | Yes | **Best for double** (Vitis DSP Lib) |
| `hls::x_complex<>` | N/A | N/A | **Deprecated** (Vivado HLS era, removed) |

### Option 1: Fixed-Point (Preferred for FPGA efficiency)
```cpp
#include <complex>
#include "ap_fixed.h"

// Limit to 27 bits total (int + frac) for optimal DSP48 block mapping
typedef ap_fixed<27, 8> fixed_t;
typedef std::complex<fixed_t> complex_t;
```

### Option 2: Floating-Point with Vitis DSP Library
```cpp
// Requires: clone Vitis_Libraries, add to include path:
// REPO_PATH/dsp/L1/include/hw/vitis_2dfft/float/
#include "vt_fft.hpp"

using namespace xf::dsp::fft;
typedef float T_innerData;
typedef complex_wrapper<T_innerData> T_elemType;
```

### Option 3: std::complex<float> with Manual Wrappers
If you must use `std::complex<float>`, wrap all operations in static functions with `INLINE off`:

```cpp
#include <complex>
using data_t    = float;
using complex_t = std::complex<float>;

// Accessors
static inline data_t c_re(const complex_t &z) { return z.real(); }
static inline data_t c_im(const complex_t &z) { return z.imag(); }

// Operations - INLINE off prevents template replication issues
static complex_t c_add(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)+c_re(b), c_im(a)+c_im(b));
}
static complex_t c_sub(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)-c_re(b), c_im(a)-c_im(b));
}
static complex_t c_mul(const complex_t &a, const complex_t &b) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)*c_re(b) - c_im(a)*c_im(b),
                   c_re(a)*c_im(b) + c_im(a)*c_re(b));
}
static complex_t c_scale(const complex_t &a, data_t s) {
  #pragma HLS INLINE off
  return complex_t(c_re(a)*s, c_im(a)*s);
}
static data_t c_abs2(const complex_t &a) {
  #pragma HLS INLINE off
  return c_re(a)*c_re(a) + c_im(a)*c_im(a);
}
static complex_t c_conj(const complex_t &a) {
  #pragma HLS INLINE off
  return complex_t(c_re(a), -c_im(a));
}
```

### Host/Testbench Compilation Compatibility
```cpp
// Provide hls::exp (and other hls:: math) for host compilation
namespace hls { using std::exp; using std::sqrt; using std::sin; using std::cos; }
```

---

## HLS Pragmas Reference

### Interface Pragmas (Top Function)
```cpp
// AXI-Lite control interface
#pragma HLS INTERFACE s_axilite port=return bundle=control

// Memory-mapped AXI for arrays (DDR access)
#pragma HLS INTERFACE m_axi port=data_in offset=slave bundle=gmem depth=1024
#pragma HLS INTERFACE s_axilite port=data_in bundle=control

// Scalar parameters via AXI-Lite
#pragma HLS INTERFACE s_axilite port=param bundle=control

// Streaming interface
#pragma HLS INTERFACE axis port=stream_in
#pragma HLS INTERFACE axis port=stream_out
```

### Storage Pragmas
```cpp
// Dual-port BRAM (default for large arrays)
#pragma HLS BIND_STORAGE variable=arr type=ram_2p impl=bram

// Single-port BRAM
#pragma HLS BIND_STORAGE variable=arr type=ram_1p impl=bram

// LUTRAM (small arrays, faster but uses LUTs)
#pragma HLS BIND_STORAGE variable=arr type=ram_1p impl=lutram

// Registers (very small arrays)
#pragma HLS BIND_STORAGE variable=arr type=ram_1p impl=auto
```

### Function Pragmas
```cpp
#pragma HLS INLINE off    // Prevent inlining (share hardware resources)
#pragma HLS INLINE        // Force inlining (duplicate for parallelism)
#pragma HLS PIPELINE II=1 // Pipeline function, initiation interval = 1
#pragma HLS DATAFLOW      // Enable task-level parallelism between functions
```

### Loop Pragmas
```cpp
#pragma HLS PIPELINE II=1              // Pipeline loop iterations
#pragma HLS UNROLL                     // Fully unroll loop
#pragma HLS UNROLL factor=4            // Partially unroll
#pragma HLS LOOP_TRIPCOUNT min=1 max=100  // Hint for variable bounds
#pragma HLS LOOP_FLATTEN               // Flatten nested loops
```

### Array Pragmas
```cpp
#pragma HLS ARRAY_PARTITION variable=arr complete dim=1   // Full partition
#pragma HLS ARRAY_PARTITION variable=arr cyclic factor=4 dim=1
#pragma HLS ARRAY_PARTITION variable=arr block factor=4 dim=1
#pragma HLS ARRAY_RESHAPE variable=arr complete dim=1     // Reshape for wider access
```

---

## Arbitrary Precision Types

```cpp
#include "ap_int.h"    // Arbitrary precision integers
#include "ap_fixed.h"  // Arbitrary precision fixed-point

// Integer types
ap_int<17> signed_17bit;
ap_uint<24> unsigned_24bit;

// Fixed-point: ap_fixed<TOTAL_BITS, INTEGER_BITS>
ap_fixed<16, 4> signed_fixed;    // 16 bits total, 4 integer, 12 fractional
ap_ufixed<16, 4> unsigned_fixed;

// For DSP48 optimization: limit to 27 bits total
ap_fixed<27, 8> dsp_optimized;  // Maps to single DSP block multiply
```

---

## Fixed-Point Implementation Best Practices

Converting floating-point algorithms to fixed-point requires careful consideration of range, precision, and HLS math function limitations.

### Type Selection Strategy

```cpp
#include "ap_fixed.h"
#include "hls_math.h"

// Tunable parameters - adjust based on your value ranges
#define DATA_TOTAL_BITS  27   // Total bits (int + frac)
#define DATA_INT_BITS    8    // Integer bits (determines range)

// Main arithmetic type - AP_RND rounds, AP_SAT saturates on overflow
typedef ap_fixed<DATA_TOTAL_BITS, DATA_INT_BITS, AP_RND, AP_SAT> data_t;

// Range for ap_fixed<27,8>:
//   Max: +127.9999998 (2^7 - 2^-19)
//   Min: -128
//   Resolution: 2^-19 ≈ 1.9e-6
```

### Quantization Mode Reference

| Mode | Behavior | Use Case |
|------|----------|----------|
| `AP_RND` | Round to nearest | General arithmetic (default choice) |
| `AP_TRN` | Truncate toward -∞ | Faster, less accurate |
| `AP_SAT` | Saturate on overflow | Prevents wraparound errors |
| `AP_WRAP` | Wrap on overflow | Faster but can cause large errors |

**Recommended combination:** `ap_fixed<W, I, AP_RND, AP_SAT>` for safety.

### HLS Math Function Bit-Width Limitations (CRITICAL)

**Problem:** `hls::exp()`, `hls::sin()`, `hls::cos()`, `hls::sqrt()` only support specific bit-widths.

| Function | Supported Types |
|----------|----------------|
| `hls::exp()` | `ap_fixed<16,8>` or narrower |
| `hls::sin()`, `hls::cos()` | `ap_fixed<32,I>` where I ≤ 3 |
| `hls::sqrt()` | Up to `ap_fixed<49,I>` |

**Solution:** Create wrapper functions that narrow input, call HLS function, then widen output:

```cpp
// Narrow type for hls::exp() compatibility
typedef ap_fixed<16, 8, AP_RND, AP_SAT> exp_t;

// Wrapper: narrow → compute → widen
static data_t fp_exp(data_t x) {
    #pragma HLS INLINE off
    exp_t x_narrow = (exp_t)x;      // Narrow to 16,8
    exp_t result = hls::exp(x_narrow);
    return (data_t)result;          // Widen back to 27,8
}
```

### Precomputing Constants That Overflow

**Problem:** Physical constants often exceed fixed-point range.
Example: wavenumber k = 7.85×10⁶ overflows `ap_fixed<27,8>` (max ≈ 128).

**Solution:** Precompute combined coefficients that fit in range:

```cpp
// BAD - k overflows ap_fixed<27,8>:
// static const data_t k = 7853981.6339f;
// data_t coeff = dz / (4 * k * dx * dx);  // Runtime overflow!

// GOOD - precompute the combined coefficient:
// dz/(4·k·dx²) = 1e-4 / (4 × 7.85e6 × (7.03e-7)²) ≈ 6.44
static const data_t UNG_COEFF = data_t(6.4368);

// Use precomputed value directly
complex_t ung = c_make(data_t(0), UNG_COEFF);
```

### Fixed-Point Testbench Considerations

1. **Relaxed error thresholds:** Use 1e-2 instead of 1e-3 for RMS error
2. **Range analysis:** Log min/max/mean of inputs and outputs
3. **Overflow detection:** Check for values near saturation limits
4. **Effective bits metric:** Compare `−log₂(relative_error)` to expected precision

```cpp
// Detect potential saturation
double sat_max = (1 << (DATA_INT_BITS - 1)) - 1;  // e.g., 127 for 8 int bits
if (std::abs(value) > sat_max * 0.99) {
    overflow_count++;
}

// Effective bits of accuracy
double eff_bits = -std::log2(std::abs((computed - golden) / golden));
// For ap_fixed<27,8>: expect ~19 bits (fractional bits)
```

### Migration Checklist (Float → Fixed)

1. **Analyze value ranges** in your algorithm (min, max for all variables)
2. **Choose bit-widths:** Integer bits ≥ ceil(log₂(max_value)) + 1 (sign)
3. **Identify HLS math calls** and check supported bit-widths
4. **Create wrapper functions** for unsupported bit-widths
5. **Precompute large constants** as combined coefficients
6. **Replace float literals** with `data_t(value)` constructors
7. **Update testbench** with range logging and relaxed thresholds
8. **Run C simulation** to verify numerical accuracy before synthesis

### References - Fixed-Point in Vitis HLS

- [UG1399: Arbitrary Precision Types](https://docs.amd.com/r/en-US/ug1399-vitis-hls/Using-Arbitrary-Precision-Data-Types)
- [UG1399: HLS Math Library](https://docs.amd.com/r/en-US/ug1399-vitis-hls/HLS-Math-Library)
- [UG1399: Fixed-Point Math Functions](https://docs.amd.com/r/en-US/ug1399-vitis-hls/Fixed-Point-Math-Functions)
- [Xilinx AR#69908: hls::exp supported types](https://support.xilinx.com/s/article/69908)

---

## Common Synthesis Issues & Solutions

### LUT Overflow (>100%)
1. Use BRAM instead of LUTRAM: `#pragma HLS BIND_STORAGE variable=x impl=bram`
2. Add `#pragma HLS INLINE off` to prevent function replication
3. Reduce `UNROLL` factors
4. Increase pipeline II (trade throughput for area)

### Long Synthesis Time
1. Add `#pragma HLS LOOP_TRIPCOUNT` for variable-bound loops
2. Use `#pragma HLS ALLOCATION` to limit resource instances
3. Simplify deeply nested loops

### Timing Failures
1. Reduce clock frequency target
2. Add pipeline stages
3. Use `#pragma HLS LATENCY min=X max=Y` constraints

### II (Initiation Interval) Not Met
1. Check for loop-carried dependencies
2. Partition arrays accessed in pipeline
3. Use `#pragma HLS DEPENDENCE` to clarify false dependencies

---

## Project File Structure (Convention)

```
ComponentName/
├── component_name.cpp      # Main HLS source (synthesized)
├── component_name.h        # Header: types, constants, top function declaration
├── component_name_tb.cpp   # Testbench (host simulation, not synthesized)
├── hls_config.cfg          # Vitis HLS project configuration (optional)
└── compile_commands.json   # IDE support (optional)
```

---

## Testbench Guidelines

Testbenches (`*_tb.cpp`) run on the host CPU and CAN use full C++:
- `printf`, `std::cout` for debugging
- STL containers for golden reference
- File I/O for test vectors
- Compare HLS output against reference with tolerance

```cpp
// Typical testbench structure
int main() {
    // 1. Initialize test data
    // 2. Call DUT (Device Under Test) - your top function
    // 3. Compare against golden reference
    // 4. Return 0 for pass, non-zero for fail
    return (error_count == 0) ? 0 : 1;
}
```

---

## References

- [Vitis HLS User Guide (UG1399)](https://docs.amd.com/r/en-US/ug1399-vitis-hls)
- [Vitis HLS Pragmas Reference](https://docs.amd.com/r/en-US/ug1399-vitis-hls/HLS-Pragmas)
- [Vitis DSP Library - GitHub](https://github.com/Xilinx/Vitis_Libraries/tree/master/dsp)
- [Vitis HLS Migration Guide](https://docs.amd.com/r/en-US/ug1399-vitis-hls/Vitis-HLS-Migration-Guide)
- [Deprecated Features](https://docs.amd.com/r/en-US/ug1399-vitis-hls/Deprecated-and-Unsupported-Features)
- [Vitis Tutorials - GitHub](https://github.com/Xilinx/Vitis-Tutorials)

---
---
# PART 2: PROJECT-SPECIFIC CONFIGURATION
---
---

## This Workspace: BPM Laser Simulation Accelerators

**Target Platform:** Kria KV260 (Zynq UltraScale+ MPSoC)

### Physical Application
BPM (Beam Propagation Method) simulation for laser beam propagation in nonlinear media using the **Crank-Nicolson ADI split-step scheme**.

**Physical effects implemented:**
- Linear diffraction (paraxial wave equation)
- Linear absorption (alpha parameter, units: m⁻¹)
- Two-photon absorption / TPA (beta parameter, units: m/W)

**Grid:** 64x64 complex field (`complex_t phi[N][N]` where `N=64`)

---

## Component Hierarchy

| Component | Description | Status |
|-----------|-------------|--------|
| **Diff_Losses_Nsteps** | Full BPM with Strang splitting for absorption | Working, synthesized |
| **Diff_Abs_Nsteps** | Full BPM with simpler absorption scheme | 101% LUT - needs optimization |
| **Diff_Abs_Nsteps_opt1** | Optimization attempt of above | In progress |
| **DiffractionOnly_Nsteps** | Diffraction only (no absorption) | Working |
| **DiffractionOnly_v1** | Earlier diffraction version | Legacy |
| **Custom_Thomas_Solver_v1** | Standalone Thomas solver | Utility |
| **ComplexMultiplier** | Complex multiplication IP | Utility |
| **BufferOverlay** | Simple buffer passthrough | Baseline test |

---

## Key Algorithms in This Project

### Thomas Solver (Tridiagonal)
Solves `Ax = b` where A is tridiagonal matrix. Used for ADI implicit steps.
- Forward elimination + backward substitution
- Purely iterative (no recursion)
- Handles boundary conditions via modified diagonal entries

### ADI Split-Step Flow
```
for s = 0 to steps-1:
    adi_x(phi, tmp)              // Implicit half-step in X direction
    apply_absorption(tmp)         // Nonlinear terms at midpoint
    adi_y(tmp, phi)              // Implicit half-step in Y direction
```

### Strang Splitting (in Diff_Losses_Nsteps)
Symmetric application for higher-order accuracy:
```
apply_N_midpoint(phi, alpha, beta):
    half_linear_absorption(phi, alpha)     // exp(-alpha * dz/4)
    half_2photon_absorption(phi, beta)     // exp(-beta * dz/4 * |phi|²)
    half_2photon_absorption(phi, beta)     // exp(-beta * dz/4 * |phi|²)
    half_linear_absorption(phi, alpha)     // exp(-alpha * dz/4)
```

---

## Physical Constants (in code)

```cpp
static constexpr data_t Lx   = 45e-6f;      // Domain size X (45 μm)
static constexpr data_t Ly   = 45e-6f;      // Domain size Y (45 μm)
static constexpr data_t dx   = Lx / N;      // Grid spacing X
static constexpr data_t dy   = Ly / N;      // Grid spacing Y
static constexpr data_t dz   = 1e-4f;       // Propagation step (100 μm)
static constexpr data_t k    = 7853981.6339f;  // Wavenumber (~2π/0.8μm for 800nm)
```

---

## Top Function Signature (Diff_Losses_Nsteps)

```cpp
void diff_losses(
    const complex_t phi_in[N][N],   // Input field (m_axi gmem)
          complex_t phi_out[N][N],  // Output field (m_axi gmem)
    int steps,                       // Number of propagation steps (s_axilite)
    data_t alpha,                    // Linear absorption coeff (s_axilite)
    data_t beta                      // TPA coefficient (s_axilite)
);
```

---

## Resource Notes for Kria KV260

- Recent synthesis of Diff_Abs_Nsteps hit **101% LUT** - not implementable
- Solution: Adjust pragmas, use more BRAM, reduce parallelism
- Internal matrices use `#pragma HLS BIND_STORAGE type=ram_2p impl=bram`
- Complex wrapper functions use `#pragma HLS INLINE off` to prevent replication
