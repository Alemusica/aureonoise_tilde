#pragma once

#include "aureo_math.hpp"

namespace aureo {

struct FieldState {
  struct EnvelopeShape {
    double attackEnd = 0.12;
    double decayEnd = 0.32;
    double releaseStart = 0.82;
    double sustainLevel = 0.55;
  };

  double temperature = 0.45;
  double last_gap_samples = 0.0;
  double last_dur_samples = 0.0;
  double last_center_distance = 0.0;
  double last_io_ratio = 1.0;

  EnvelopeShape make_envelope(double attack_ratio,
                              double decay_ratio,
                              double sustain_level,
                              double release_ratio,
                              double gap_samples,
                              double dur_samples,
                              double center_distance)
  {
    sustain_level = clamp01(sustain_level);
    attack_ratio = clamp(attack_ratio, 0.01, 4.0);
    decay_ratio = clamp(decay_ratio, 0.01, 4.0);
    release_ratio = clamp(release_ratio, 0.01, 4.0);
    gap_samples = std::max(0.0, gap_samples);
    dur_samples = std::max(1.0, dur_samples);
    center_distance = clamp01(center_distance);

    last_gap_samples = gap_samples;
    last_dur_samples = dur_samples;
    last_center_distance = center_distance;

    const double io_ratio = gap_samples / dur_samples;
    last_io_ratio = io_ratio;

    const double io_weight = clamp(0.15 + 0.7 * (io_ratio / (io_ratio + 1.0)), 0.15, 0.95);
    const double center_boost = 0.55 + 0.45 * (1.0 - center_distance);
    const double dynamic_span = clamp(io_weight * center_boost, 0.15, 0.95);

    const double total = attack_ratio + decay_ratio + release_ratio;
    const double safe_total = (total <= 1.0e-9) ? 1.0 : total;
    double attackNorm = (attack_ratio / safe_total) * dynamic_span;
    double decayNorm = (decay_ratio / safe_total) * dynamic_span;
    double releaseNorm = (release_ratio / safe_total) * dynamic_span;

    double sustainNorm = 1.0 - (attackNorm + decayNorm + releaseNorm);
    if (sustainNorm < 0.05) {
      const double deficit = 0.05 - sustainNorm;
      const double sumADR = attackNorm + decayNorm + releaseNorm;
      if (sumADR > 1.0e-6) {
        const double scale = 1.0 - deficit / sumADR;
        attackNorm *= scale;
        decayNorm *= scale;
        releaseNorm *= scale;
      }
      sustainNorm = 0.05;
    }

    EnvelopeShape env;
    env.attackEnd = clamp(attackNorm, 1.0e-4, 0.75);
    env.decayEnd = clamp(env.attackEnd + std::max(1.0e-4, decayNorm), env.attackEnd + 1.0e-4, 0.96);
    const double relStart = 1.0 - clamp(releaseNorm, 1.0e-4, 0.95);
    env.releaseStart = clamp(std::max(env.decayEnd, relStart), env.decayEnd + 1.0e-4, 0.999);
    env.sustainLevel = sustain_level;
    return env;
  }

  double envelope(double phase, const EnvelopeShape& env) const
  {
    phase = clamp01(phase);
    if (phase <= env.attackEnd) {
      const double denom = std::max(1.0e-6, env.attackEnd);
      return phase / denom;
    }
    if (phase <= env.decayEnd) {
      const double t = (phase - env.attackEnd) / std::max(1.0e-6, env.decayEnd - env.attackEnd);
      return 1.0 + (env.sustainLevel - 1.0) * t;
    }
    if (phase < env.releaseStart) {
      return env.sustainLevel;
    }
    const double t = (phase - env.releaseStart) / std::max(1.0e-6, 1.0 - env.releaseStart);
    return env.sustainLevel * (1.0 - t);
  }
};

} // namespace aureo

