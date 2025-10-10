#pragma once

#include "aureo_math.hpp"

namespace aureo {

struct BinauralCoefficients {
  double gainL = 1.0;
  double gainR = 1.0;
  double itd_samples = 0.0;
  double ild_db = 0.0;
  double azimuth_rad = 0.0;
  double focus = 1.0;
  double crossfeed = 0.0;
};

inline void pan_equal_power(double pan, double width, double& gL, double& gR)
{
  const double p = clamp(pan * width, -1.0, 1.0);
  const double th = (p + 1.0) * (kPi * 0.25);
  gL = std::cos(th);
  gR = std::sin(th);
}

inline double woodworth_itd(double theta)
{
  const double limit = 0.5 * kPi;
  if (theta > limit) theta = limit;
  if (theta < -limit) theta = -limit;
  return theta + std::sin(theta);
}

inline double map_itd_samples_frac(double sr, double itd_us, double pan, double u)
{
  sr = (sr > 0.0) ? sr : 44100.0;
  const double maxSec = clamp(itd_us, 0.0, 1600.0) * 1.0e-6;
  if (maxSec <= 0.0) return 0.0;

  const double headRadius = 0.0875; // ~17.5 cm diameter
  const double speedOfSound = 343.0;
  const double theta = clamp(pan, -1.0, 1.0) * (kPi * 0.5);
  const double geom = woodworth_itd(theta);
  const double geomMax = woodworth_itd(0.5 * kPi);
  double itdSec = (headRadius / speedOfSound) * geom;
  const double itdSecMax = (headRadius / speedOfSound) * geomMax;
  if (itdSecMax > 1.0e-9) {
    itdSec = (maxSec / itdSecMax) * itdSec;
  } else {
    itdSec = 0.0;
  }

  const double jitter = (2.0 * clamp01(u) - 1.0) * 0.08;
  itdSec = clamp(itdSec + jitter * maxSec, -maxSec, maxSec);
  return itdSec * sr;
}

inline double map_ild_db(double ild_db_max, double pan, double u)
{
  const double maxDb = clamp(ild_db_max, 0.0, 30.0);
  if (maxDb <= 1.0e-9) return 0.0;

  const double theta = clamp(pan, -1.0, 1.0) * (kPi * 0.5);
  const double lateral = std::sin(theta);
  const double lateralMag = std::pow(std::abs(lateral), 0.9);
  const double random = map_phi_range(0.35, 1.0, clamp01(std::abs(2.0 * u - 1.0)));
  double ild = maxDb * lateralMag * random;
  const double frontFocus = 0.65 + 0.35 * lateralMag;
  ild *= frontFocus;
  return (lateral >= 0.0) ? ild : -ild;
}

inline BinauralCoefficients compute_binaural(double sr,
                                             double pan,
                                             double width,
                                             double itd_us,
                                             double ild_db_max,
                                             double itd_rand,
                                             double ild_rand)
{
  BinauralCoefficients coeffs;
  const double clampedWidth = clamp(width, 0.0, 1.0);
  const double panNorm = clamp(pan * clampedWidth, -1.0, 1.0);

  pan_equal_power(pan, clampedWidth, coeffs.gainL, coeffs.gainR);
  coeffs.azimuth_rad = panNorm * (kPi * 0.5);

  const double lateral = std::sin(coeffs.azimuth_rad);
  coeffs.focus = 0.25 + 0.75 * std::pow(std::abs(lateral), 1.35);
  coeffs.crossfeed = (1.0 - coeffs.focus) * 0.18;

  coeffs.itd_samples = map_itd_samples_frac(sr, itd_us, panNorm, itd_rand);
  coeffs.ild_db = map_ild_db(ild_db_max, panNorm, ild_rand);
  return coeffs;
}

} // namespace aureo

