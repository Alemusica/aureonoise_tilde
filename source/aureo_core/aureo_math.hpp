#pragma once

#include <algorithm>
#include <cmath>

#include "aureo_config.hpp"

namespace aureo {

inline double clamp(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
inline double clamp01(double x) { return clamp(x, 0.0, 1.0); }
inline double db_to_lin(double dB) { return std::pow(10.0, dB / 20.0); }
inline double hann01(double p) { p = clamp01(p); return 0.5 - 0.5 * std::cos(2.0 * kPi * p); }
inline double soft_tanh(double v) { return std::tanh(v); }

inline double map_phi_range(double vmin, double vmax, double u)
{
  u = clamp01(u);
  vmin = std::max(1.0e-12, vmin);
  vmax = std::max(vmin * 1.000001, vmax);
  const double K = std::log(vmax / vmin) / std::log(kPhi);
  return vmin * std::pow(kPhi, K * u);
}

inline int map_sr_hold_base(double amt, double u)
{
  amt = clamp01(amt);
  double steps = map_phi_range(1.0, 64.0, amt);
  int N = (int)std::max(1.0, std::round(steps));
  int jitter = (int)std::floor(clamp01(u) * 0.999 * N);
  return std::max(1, N - jitter);
}

inline int map_bits(double amt)
{
  amt = clamp01(amt);
  int bits = (int)std::round(16.0 - amt * 12.0);
  return std::max(4, std::min(16, bits));
}

inline double quantize_bits(double x, int bits)
{
  int levels = (1 << (bits - 1)) - 1;
  if (levels <= 0) return 0.0;
  return std::round(x * levels) / (double)levels;
}

} // namespace aureo

