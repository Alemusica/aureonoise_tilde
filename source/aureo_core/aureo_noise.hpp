#pragma once

#include <cstdint>

#include "aureo_math.hpp"
#include "aureo_rng.hpp"

namespace aureo {

enum class NoiseColor : int { White = 0, Pink, Brown };

enum class GrainKind : int { Burst = 0, VhsDrop, Stutter, Aliaser };

struct NoiseColorState {
  NoiseColor color = NoiseColor::Pink;
  double amount = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;

  void reset() { z1 = 0.0; z2 = 0.0; }

  double process(RNG& rng)
  {
    const double w = rng.uniPM1();
    const double amt = clamp01(amount);
    if (color == NoiseColor::White || amt <= 1.0e-6) return w;

    if (color == NoiseColor::Pink) {
      z1 = (1.0 - 0.02 * amt) * z1 + (0.02 * amt) * w;
      return (1.0 - amt) * w + amt * z1;
    } else {
      z2 = (1.0 - 0.001 * amt) * z2 + (0.03 * amt) * w;
      return clamp((1.0 - amt) * w + amt * z2, -1.5, 1.5);
    }
  }
};

struct Biquad {
  double b0 = 1.0, b1 = 0.0, b2 = 0.0;
  double a1 = 0.0, a2 = 0.0;
  double z1 = 0.0, z2 = 0.0;

  double process(double x)
  {
    double y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
  }

  void reset() { z1 = z2 = 0.0; }
};

struct Grain {
  bool on = false;
  uint32_t age = 0;
  uint32_t dur = 0;
  double amp = 0.0;
  double panL = 1.0;
  double panR = 1.0;
  double itd = 0.0;
  double gL = 1.0;
  double gR = 1.0;
  int sr_holdN = 1;
  int sr_holdCnt = 1;
  double heldL = 0.0;
  double heldR = 0.0;
  int q_levels = 0;
  GrainKind kind = GrainKind::Burst;
};

inline GrainKind choose_kind(double mix, double u)
{
  const double m = clamp01(mix);
  if (u < 0.25 * m) return GrainKind::VhsDrop;
  if (u < 0.60 * m) return GrainKind::Stutter;
  if (u < 1.00 * m) return GrainKind::Aliaser;
  return GrainKind::Burst;
}

} // namespace aureo

