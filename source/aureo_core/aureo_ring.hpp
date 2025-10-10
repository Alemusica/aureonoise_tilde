#pragma once

#include <array>
#include <cstdint>
#include <cmath>

#include "aureo_config.hpp"

namespace aureo {

struct RingState {
  std::array<double, kRingSize> data{};
  uint32_t writeIndex = 0u;

  void clear()
  {
    data.fill(0.0);
    writeIndex = 0u;
  }
};

inline double lagrange3(const double* r, double pos)
{
  const int64_t i = (int64_t)std::floor(pos);
  const double  f = pos - (double)i;
  const uint32_t i_1 = (uint32_t)(i - 1) & kRingMask;
  const uint32_t i0  = (uint32_t)(i + 0) & kRingMask;
  const uint32_t i1  = (uint32_t)(i + 1) & kRingMask;
  const uint32_t i2  = (uint32_t)(i + 2) & kRingMask;

  const double x_1 = r[i_1], x0 = r[i0], x1 = r[i1], x2 = r[i2];
  const double f1 = f - 1.0;
  const double f2 = f - 2.0;

  const double c_1 = -f * (f1) * (f2) / 6.0;
  const double c0  = (f + 1.0) * (f1) * (f2) / 2.0;
  const double c1  = -f * (f + 1.0) * (f2) / 2.0;
  const double c2  = f * (f + 1.0) * f1 / 6.0;

  return c_1 * x_1 + c0 * x0 + c1 * x1 + c2 * x2;
}

inline void ring_read_stereo_itd_frac(const RingState& ring, uint32_t base, double d, double& oL, double& oR)
{
  const double ad = std::abs(d);
  const double B  = std::ceil(ad) + 2.0;
  const double DL = B + (d > 0.0 ? d : 0.0);
  const double DR = B + (d < 0.0 ? -d : 0.0);
  const double posL = (double)base - DL;
  const double posR = (double)base - DR;
  oL = lagrange3(ring.data.data(), posL);
  oR = lagrange3(ring.data.data(), posR);
}

} // namespace aureo

