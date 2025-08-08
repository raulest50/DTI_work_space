// hls_stub.h
#ifndef HLS_STUB_H
#define HLS_STUB_H
#include <complex>
#include <cmath>
#include <queue>
namespace hls {
template<typename T>
using x_complex = std::complex<T>;

// ---------------------------------------------------------------------
// Minimal hls::stream stub used for host compilation
// ---------------------------------------------------------------------
template<typename T>
class stream {
    std::queue<T> q;
public:
    inline bool empty() const { return q.empty(); }
    inline void write(const T& v) { q.push(v); }
    inline T read() { T v = q.front(); q.pop(); return v; }
};

inline float exp(float x) { return std::exp(x); }
inline double exp(double x) { return std::exp(x); }

inline float hypot(float x, float y) { return std::hypot(x, y); }
inline double hypot(double x, double y) { return std::hypot(x, y); }

inline void sincos(float x, float* s, float* c) { *s = std::sin(x); *c = std::cos(x); }
inline void sincos(double x, double* s, double* c) { *s = std::sin(x); *c = std::cos(x); }
}
#endif
