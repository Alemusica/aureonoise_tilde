#pragma once

#include "aureo_math.hpp"

namespace aureo {

inline void pan_equal_power(double pan, double width, double& gL, double& gR)
{
  const double p = clamp(pan * width, -1.0, 1.0);
  const double th = (p + 1.0) * (kPi * 0.25);
  gL = std::cos(th);
  gR = std::sin(th);
}

inline double map_itd_samples_frac(double sr, double itd_us, double u)
{
  const double maxSamp = clamp(itd_us, 0.0, 1000.0) * 1.0e-6 * sr;
#if AUREO_ITD_PHI_SHAPE
  const double m    = std::abs(2.0 * u - 1.0);
  const double mag  = map_phi_range(0.0, std::max(0.0, maxSamp), m);
  const double sign = (u < 0.5) ? -1.0 : +1.0;
  return sign * mag;
#else
  return (2.0 * u - 1.0) * maxSamp;
#endif
}

inline double map_ild_gain(double itd_db_max, double u, bool left)
{
  const double Dmax = clamp(itd_db_max, 0.0, 24.0);
  const double m    = std::abs(2.0 * u - 1.0);
  const double dBmag= map_phi_range(0.0, Dmax, m);
  const double sgn  = (u < 0.5 ? +1.0 : -1.0);
  const double dL   = left ? (sgn * dBmag) : (-sgn * dBmag);
  return db_to_lin(dL);
}

} // namespace aureo

