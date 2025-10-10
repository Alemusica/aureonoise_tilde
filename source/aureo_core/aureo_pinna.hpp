#pragma once

#include <utility>

#include "aureo_binaural.hpp"

namespace aureo {

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

  double itd_samples(double sr, double u) const { return map_itd_samples_frac(sr, itd_us, u); }

  double ild_gain(bool left, double u) const { return map_ild_gain(ild_db, u, left); }
};

} // namespace aureo

