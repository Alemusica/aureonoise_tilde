#pragma once

#include <complex>
#include <cstdint>

#include "aureo_math.hpp"
#include "aureo_rng.hpp"

namespace aureo {

struct FieldState {
  struct EnvelopeShape {
    double attackEnd = 0.12;
    double decayEnd = 0.32;
    double releaseStart = 0.82;
    double sustainLevel = 0.55;
  };

  struct ComplexEnvelopeParams {
    double magnitudeVariance = 0.18;
    double phaseVariance = 0.25;
    double correlation = 0.72;
    double tauMs = 24.0;
    double sustainBias = 0.35;
  };

  struct ComplexEnvelopeState {
    std::complex<double> value{0.0, 0.0};
    double magNoise = 0.0;
    double phaseNoise = 0.0;
    double response = 0.05;
    double magnitudeVariance = 0.0;
    double phaseVariance = 0.0;
    double correlation = 0.0;
    double sustainBias = 0.0;
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

  ComplexEnvelopeParams make_complex_params(double magnitudeVariance,
                                            double phaseVariance,
                                            double correlation,
                                            double tauMs,
                                            double sustainBias) const
  {
    ComplexEnvelopeParams params;
    params.magnitudeVariance = clamp(magnitudeVariance, 0.0, 2.0);
    params.phaseVariance = clamp(phaseVariance, 0.0, 2.0);
    params.correlation = clamp(correlation, 0.0, 0.9995);
    params.tauMs = clamp(tauMs, 0.1, 4000.0);
    params.sustainBias = clamp01(sustainBias);
    return params;
  }

  ComplexEnvelopeState prepare_complex_envelope(uint32_t dur_samples,
                                                double sr,
                                                RNG& rng,
                                                const EnvelopeShape& env,
                                                const ComplexEnvelopeParams& params) const
  {
    ComplexEnvelopeState state;
    state.magnitudeVariance = std::max(0.0, params.magnitudeVariance);
    state.phaseVariance = std::max(0.0, params.phaseVariance);
    state.correlation = clamp(params.correlation, 0.0, 0.9995);

    const double sustain = clamp01(env.sustainLevel);
    state.sustainBias = clamp01(params.sustainBias) * (0.25 + 0.75 * sustain);

    const double sr_safe = std::max(1.0, sr);
    double tau_samples = params.tauMs * 0.001 * sr_safe;
    tau_samples *= 0.6 + 0.4 * sustain;
    if (!std::isfinite(tau_samples) || tau_samples <= 0.5) {
      tau_samples = 0.25 * std::max(8.0, static_cast<double>(dur_samples));
    }
    const double denom = std::max(1.0, tau_samples);
    state.response = 1.0 - std::exp(-1.0 / denom);

    const double initPhase = 2.0 * kPi * rng.uni01();
    state.value = std::polar(0.0, initPhase);
    state.magNoise = 0.0;
    state.phaseNoise = 0.0;

    return state;
  }

  double advance_complex_envelope(double phase,
                                  const EnvelopeShape& env,
                                  ComplexEnvelopeState& state,
                                  RNG& rng) const
  {
    phase = clamp01(phase);
    const double base = envelope(phase, env);

    const double corr = clamp(state.correlation, 0.0, 0.9995);
    const double one_minus_corr = 1.0 - corr;
    const double magNoise = corr * state.magNoise + one_minus_corr * rng.uniPM1();
    const double phaseNoise = corr * state.phaseNoise + one_minus_corr * rng.uniPM1();
    state.magNoise = magNoise;
    state.phaseNoise = phaseNoise;

    double target = base * (1.0 + state.magnitudeVariance * magNoise);
    if (state.sustainBias > 1.0e-6) {
      target = (1.0 - state.sustainBias) * target + state.sustainBias * env.sustainLevel;
    }
    target = clamp(target, 0.0, 1.5);

    const double currentMag = std::abs(state.value);
    const double slew = clamp(state.response, 1.0e-5, 1.0);
    const double newMag = currentMag + (target - currentMag) * slew;

    double newPhase = std::arg(state.value);
    if (!std::isfinite(newPhase)) newPhase = 0.0;
    const double phaseScale = state.phaseVariance * (0.25 + 0.75 * base);
    newPhase += phaseScale * phaseNoise;

    state.value = std::polar(std::max(0.0, newMag), newPhase);
    return std::abs(state.value);
  }
};

} // namespace aureo

