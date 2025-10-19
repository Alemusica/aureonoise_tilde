#pragma once

#include <cstdint>

#include "aureo_math.hpp"
#include "aureo_rng.hpp"
#include "aureo_field_adsr.hpp"

namespace aureo {

enum class NoiseColor : int { White = 0, Pink, Brown };

enum class GrainKind : int { Burst = 0, VhsDrop, Stutter, Aliaser };

struct NoiseColorState {
  NoiseColor color = NoiseColor::Pink;
  double amount = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
  double z3 = 0.0;

  void reset() { z1 = 0.0; z2 = 0.0; z3 = 0.0; }

  double process(RNG& rng)
  {
    const double w = rng.uniPM1();
    const double amt = clamp01(amount);
    if (color == NoiseColor::White || amt <= 1.0e-6) return w;

    if (color == NoiseColor::Pink) {
      z1 = (1.0 - 0.02 * amt) * z1 + (0.02 * amt) * w;
      return (1.0 - amt) * w + amt * z1;
    } else {
      const double slowPole = map_phi_range(0.88, 0.995, amt);
      const double slowerPole = map_phi_range(0.70, 0.985, amt);
      z1 = slowPole * z1 + (1.0 - slowPole) * w;
      z2 = slowerPole * z2 + (1.0 - slowerPole) * z1;
      z3 = (0.5 + 0.5 * (1.0 - amt)) * z3 + (0.5 * amt) * z2;
      const double brown = z2 + 0.35 * z3;
      const double shaped = soft_tanh(brown * (1.0 + 1.4 * amt));
      return (1.0 - amt) * w + amt * shaped;
    }
  }
};

struct Grain {
  bool on = false;
  uint32_t age = 0;
  uint32_t dur = 0;
  double amp = 0.0;
  double pan = 0.0;
  double panL = 1.0;
  double panR = 1.0;
  double itd = 0.0;
  double gL = 1.0;
  double gR = 1.0;
  double crossfeed = 0.0;
  double focus = 0.0;
  double ipd_coeff = 0.0;
  double ipd_zL = 0.0;
  double ipd_zR = 0.0;
  double shadow_a = 0.0;
  double shadow_zL = 0.0;
  double shadow_zR = 0.0;
  bool shadow_left = false;
  bool shadow_right = false;
  int sr_holdN = 1;
  int sr_holdCnt = 1;
  double heldL = 0.0;
  double heldR = 0.0;
  int q_levels = 0;
  GrainKind kind = GrainKind::Burst;
  FieldState::EnvelopeShape env;
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

