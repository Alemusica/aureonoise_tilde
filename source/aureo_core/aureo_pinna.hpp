#pragma once

#include <algorithm>
#include <cmath>
#include <utility>

#include "aureo_binaural.hpp"
#include "aureo_math.hpp"

namespace aureo {

struct Biquad {
  double b0 = 1.0, b1 = 0.0, b2 = 0.0;
  double a1 = 0.0, a2 = 0.0;
  double z1 = 0.0, z2 = 0.0;

  void setNotch(double sr, double f0, double q, double minHz)
  {
    if (sr <= 0.0) sr = 44100.0;
    const double nyFactor = 0.45 * sr;
    const double clampedF0 = clamp(f0, minHz, nyFactor);
    const double normQ = std::max(0.1, q);
    const double omega = 2.0 * kPi * (clampedF0 / sr);
    const double cosw = std::cos(omega);
    const double alpha = std::sin(omega) / (2.0 * normQ);

    const double a0 = 1.0 + alpha;
    b0 = 1.0 / a0;
    b1 = (-2.0 * cosw) / a0;
    b2 = 1.0 / a0;
    a1 = (-2.0 * cosw) / a0;
    a2 = (1.0 - alpha) / a0;
  }

  double proc(double x)
  {
    const double y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
  }

  void clear() { z1 = z2 = 0.0; }
};

struct Pinna {
  double width = 1.0;
  double itd_us = 600.0;
  double ild_db = 6.0;

  std::pair<double, double> gains(double pan) const
  {
    double gL = 1.0, gR = 1.0;
    pan_equal_power(pan, width, gL, gR);
    return {gL, gR};
  }

  double itd_samples(double sr, double pan, double u) const
  {
    return map_itd_samples_frac(sr, itd_us, pan, u);
  }

  double ild_gain(bool left, double pan, double u) const
  {
    const double ild = map_ild_db(ild_db, pan, u);
    return left ? db_to_lin(ild) : db_to_lin(-ild);
  }

  BinauralCoefficients coefficients(double sr, double pan, double itd_u, double ild_u) const
  {
    return compute_binaural(sr, pan, width, itd_us, ild_db, itd_u, ild_u);
  }
};

} // namespace aureo

