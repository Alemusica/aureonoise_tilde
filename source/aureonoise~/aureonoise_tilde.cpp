// aureonoise_tilde.cpp — generatore di texture di rumore/glitch con motore generativo Aureo
// Tutto è spinto da relazioni non-periodiche: φ, √2, costante plastica, primi.
// Stereo 0in/2out, pan equal-power + ITD/ILD microbinaurale, “VHS-ish” wow/flutter, dropouts, bit/sr crush.
// Eventi granulari non ripetitivi (sequenze a bassa discrepanza tipo Weyl φ) per ridurre l’assuefazione.
// © 2025 — MIT License. Compilare C++17. Max SDK 8+. Nome oggetto in Max: aureonoise~

#include "aureonoise_state.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef __SSE2__
#include <xmmintrin.h>
#endif

static void aureonoise_reset_reports(t_aureonoise* x)
{
  if (!x) return;
  if (x->report_mu) systhread_mutex_lock(x->report_mu);
  x->report_head = 0;
  x->report_count = 0;
  x->report_total.store(0, std::memory_order_relaxed);
  x->report_dropped.store(0, std::memory_order_relaxed);
  for (auto& entry : x->report_log) entry = {};
  x->dichotic_report_head = 0;
  x->dichotic_report_count = 0;
  x->dichotic_report_total.store(0, std::memory_order_relaxed);
  x->dichotic_report_dropped.store(0, std::memory_order_relaxed);
  for (auto& entry : x->dichotic_report_log) entry = {};
  if (x->report_mu) systhread_mutex_unlock(x->report_mu);
}

static void aureonoise_log_grain_event(t_aureonoise* x, const t_aureonoise::GrainReport& report)
{
  if (!x) return;
  bool stored = false;
  if (x->report_mu) {
    if (systhread_mutex_trylock(x->report_mu) == 0) {
      const size_t cap = x->report_log.size();
      if (cap > 0) {
        x->report_log[x->report_head] = report;
        x->report_head = (x->report_head + 1) % cap;
        if (x->report_count < cap) ++x->report_count;
      }
      stored = true;
      systhread_mutex_unlock(x->report_mu);
    }
  } else {
    const size_t cap = x->report_log.size();
    if (cap > 0) {
      x->report_log[x->report_head] = report;
      x->report_head = (x->report_head + 1) % cap;
      if (x->report_count < cap) ++x->report_count;
    }
    stored = true;
  }

  if (!stored) {
    x->report_dropped.fetch_add(1, std::memory_order_relaxed);
  }
}

static void aureonoise_log_dichotic_event(t_aureonoise* x,
                                          const t_aureonoise::DichoticReport& report)
{
  if (!x) return;
  bool stored = false;
  if (x->report_mu) {
    if (systhread_mutex_trylock(x->report_mu) == 0) {
      const size_t cap = x->dichotic_report_log.size();
      if (cap > 0) {
        x->dichotic_report_log[x->dichotic_report_head] = report;
        x->dichotic_report_head = (x->dichotic_report_head + 1) % cap;
        if (x->dichotic_report_count < cap) ++x->dichotic_report_count;
      }
      stored = true;
      systhread_mutex_unlock(x->report_mu);
    }
  } else {
    const size_t cap = x->dichotic_report_log.size();
    if (cap > 0) {
      x->dichotic_report_log[x->dichotic_report_head] = report;
      x->dichotic_report_head = (x->dichotic_report_head + 1) % cap;
      if (x->dichotic_report_count < cap) ++x->dichotic_report_count;
    }
    stored = true;
  }

  if (!stored) {
    x->dichotic_report_dropped.fetch_add(1, std::memory_order_relaxed);
  }
}

static inline double aureonoise_allpass(double x, double a, double& z)
{
  const double y = -a * x + z;
  z = x + a * y;
  return y;
}

static inline double aureonoise_shadow_lp(double x, double a, double& z)
{
  const double y = (1.0 - a) * x + a * z;
  z = y;
  return y;
}

static inline double aureonoise_smoothstep01(double x)
{
  x = aureo::clamp01(x);
  return x * x * (3.0 - 2.0 * x);
}

static double aureonoise_dichotic_envelope(const t_aureonoise::DichoticRuntime& s)
{
  if (s.burst_total_samples <= 1.0) return 0.0;
  const double pos = s.burst_elapsed_samples;
  const double attack = std::max(1.0, s.attack_samples);
  const double release = std::max(1.0, s.release_samples);
  const double sustain_start = std::min(attack, s.burst_total_samples);
  double sustain_end = s.burst_total_samples - release;
  if (sustain_end < sustain_start) sustain_end = sustain_start;
  if (pos < sustain_start) {
    return aureonoise_smoothstep01(pos / sustain_start);
  }
  if (pos > sustain_end) {
    const double tail = std::max(1.0, s.burst_total_samples - sustain_end);
    return aureonoise_smoothstep01((s.burst_total_samples - pos) / tail);
  }
  return 1.0;
}

void aureonoise_reset_dichotic_state(t_aureonoise* x)
{
  if (!x) return;
  const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
  const long mode_idx = std::clamp<long>(x->p_mode, 0, 1);
  x->mode = static_cast<t_aureonoise::StimulusMode>(mode_idx);

  auto clamp_content = [](long v) {
    if (v < 0) v = 0;
    if (v > 1) v = 1;
    return v;
  };

  auto& ds = x->dichotic;
  ds.left_next = true;
  ds.active = false;
  ds.left_ear = true;
  ds.match = true;
  ds.samples_until_next = 0.0;
  ds.interval_samples = 0.0;
  ds.burst_total_samples = 0.0;
  ds.burst_elapsed_samples = 0.0;
  ds.attack_samples = 0.0;
  ds.release_samples = 0.0;
  ds.amp = aureo::clamp(x->p_dichotic_amp, 0.0, 1.0);
  ds.tone_phase = 0.0;
  ds.tone_inc = 0.0;
  ds.current_freq = 0.0;
  ds.noise_state = 0.0;
  ds.match_freq_hz = 0.0;
  ds.mismatch_freq_hz = 0.0;
  ds.match_probability = aureo::clamp01(x->p_dichotic_match_prob);
  ds.match_content = static_cast<t_aureonoise::DichoticContent>(clamp_content(x->p_dichotic_content_match));
  ds.mismatch_content = static_cast<t_aureonoise::DichoticContent>(clamp_content(x->p_dichotic_content_mismatch));
  ds.current_content = ds.match_content;
  ds.onset_timestamp_sec = 0.0;
  ds.trial_index = 0;

  const double rate = aureo::clamp(x->p_dichotic_rate_hz, 0.01, 2.0);
  ds.interval_samples = (rate > 1.0e-6) ? (sr / rate) : (sr * 10.0);

  const double burst_ms = aureo::clamp(x->p_dichotic_burst_ms, 50.0, 4000.0);
  ds.burst_total_samples = std::max(1.0, burst_ms * 0.001 * sr);

  const double min_env_ms = 5.0;
  const double max_env_ms = std::max(min_env_ms, 0.9 * burst_ms);
  const double att_ms = aureo::clamp(x->p_dichotic_env_attack_ms, min_env_ms, max_env_ms);
  double attack_ms_phi = aureo::quantize_phi(min_env_ms, max_env_ms, att_ms);
  double release_ms_phi = attack_ms_phi * aureo::kPhi;
  if (release_ms_phi > max_env_ms) {
    release_ms_phi = max_env_ms;
    attack_ms_phi = std::max(min_env_ms, release_ms_phi / aureo::kPhi);
    attack_ms_phi = aureo::quantize_phi(min_env_ms, max_env_ms, attack_ms_phi);
  }
  const double max_sum_ms = std::max(2.0 * min_env_ms, 0.98 * burst_ms);
  double env_sum_ms = attack_ms_phi + release_ms_phi;
  if (env_sum_ms > max_sum_ms) {
    const double scale = max_sum_ms / env_sum_ms;
    attack_ms_phi = std::max(min_env_ms, attack_ms_phi * scale);
    release_ms_phi = std::max(min_env_ms, release_ms_phi * scale);
    env_sum_ms = attack_ms_phi + release_ms_phi;
  }
  ds.attack_samples = std::max(1.0, attack_ms_phi * 0.001 * sr);
  ds.release_samples = std::max(1.0, release_ms_phi * 0.001 * sr);
  const double sum_env = ds.attack_samples + ds.release_samples;
  if (sum_env >= ds.burst_total_samples) {
    const double scale = (ds.burst_total_samples - 1.0) / std::max(1.0, sum_env);
    if (scale > 0.0) {
      ds.attack_samples = std::max(1.0, ds.attack_samples * scale);
      ds.release_samples = std::max(1.0, ds.release_samples * scale);
    }
  }

  const double freq_min = 50.0;
  const double freq_max = 8000.0;
  ds.match_freq_hz = aureo::quantize_phi(freq_min, freq_max, x->p_dichotic_tone_match_hz);
  ds.mismatch_freq_hz = aureo::quantize_phi(freq_min, freq_max, x->p_dichotic_tone_mismatch_hz);
  if (std::abs(ds.mismatch_freq_hz - ds.match_freq_hz) < 1.0) {
    const double candidate_up = ds.match_freq_hz * aureo::kPhi;
    const double candidate_down = ds.match_freq_hz / aureo::kPhi;
    double fallback = ds.match_freq_hz;
    if (candidate_up <= freq_max) fallback = candidate_up;
    else if (candidate_down >= freq_min) fallback = candidate_down;
    ds.mismatch_freq_hz = aureo::quantize_phi(freq_min, freq_max, fallback);
  }

  ds.samples_until_next = ds.interval_samples;
}

static void aureonoise_emit_dichotic_trial(t_aureonoise* x,
                                           const t_aureonoise::DichoticReport& report)
{
  if (!x || !x->out_info) return;
  static t_symbol* sym_dichotic = gensym("dichotic");
  static t_symbol* sym_index = gensym("index");
  static t_symbol* sym_time = gensym("time_ms");
  static t_symbol* sym_ear = gensym("ear");
  static t_symbol* sym_match = gensym("match");
  static t_symbol* sym_content = gensym("content");
  static t_symbol* sym_freq = gensym("freq_hz");
  static t_symbol* sym_amp = gensym("amp");
  static t_symbol* sym_L = gensym("L");
  static t_symbol* sym_R = gensym("R");
  static t_symbol* sym_tone = gensym("tone");
  static t_symbol* sym_noise = gensym("noise");

  t_atom av[16];
  long ac = 0;
  atom_setsym(av + ac++, sym_index);
  atom_setlong(av + ac++, static_cast<t_atom_long>(report.index));
  atom_setsym(av + ac++, sym_time);
  atom_setfloat(av + ac++, report.timestamp_sec * 1000.0);
  atom_setsym(av + ac++, sym_ear);
  atom_setsym(av + ac++, report.left_ear ? sym_L : sym_R);
  atom_setsym(av + ac++, sym_match);
  atom_setlong(av + ac++, report.match ? 1 : 0);
  atom_setsym(av + ac++, sym_content);
  atom_setsym(av + ac++, report.is_tone ? sym_tone : sym_noise);
  atom_setsym(av + ac++, sym_freq);
  atom_setfloat(av + ac++, report.freq_hz);
  atom_setsym(av + ac++, sym_amp);
  atom_setfloat(av + ac++, report.amp);
  outlet_anything(x->out_info, sym_dichotic, ac, av);
}

static void aureonoise_start_dichotic_burst(t_aureonoise* x)
{
  if (!x) return;
  auto& ds = x->dichotic;
  const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
  ds.active = true;
  ds.left_ear = ds.left_next;
  ds.left_next = !ds.left_next;
  const double u = x->rng.uni01();
  ds.match = (u <= ds.match_probability);
  ds.current_content = ds.match ? ds.match_content : ds.mismatch_content;
  const bool is_tone = (ds.current_content == t_aureonoise::DichoticContent::Tone);
  const double raw_freq = ds.match ? ds.match_freq_hz : ds.mismatch_freq_hz;
  ds.current_freq = is_tone ? aureo::clamp(raw_freq, 50.0, 8000.0) : 0.0;
  ds.tone_phase = 0.0;
  ds.tone_inc = is_tone ? (2.0 * aureo::kPi * ds.current_freq / sr) : 0.0;
  ds.noise_state = 0.0;
  ds.burst_elapsed_samples = 0.0;
  ds.onset_timestamp_sec = static_cast<double>(x->sample_counter) / sr;
  const uint64_t idx = x->dichotic_report_total.fetch_add(1, std::memory_order_relaxed) + 1;
  ds.trial_index = idx;
  t_aureonoise::DichoticReport report;
  report.index = idx;
  report.timestamp_sec = ds.onset_timestamp_sec;
  report.left_ear = ds.left_ear;
  report.match = ds.match;
  report.is_tone = is_tone;
  report.freq_hz = ds.current_freq;
  report.amp = ds.amp;
  aureonoise_log_dichotic_event(x, report);
  aureonoise_emit_dichotic_trial(x, report);
  ds.samples_until_next = ds.interval_samples;
}

static void aureonoise_perform64_dichotic(t_aureonoise* x, double* outL, double* outR, long sampleframes)
{
  auto& ds = x->dichotic;
  for (long n = 0; n < sampleframes; ++n) {
    if (!ds.active) {
      ds.samples_until_next -= 1.0;
      if (ds.samples_until_next <= 0.0) {
        aureonoise_start_dichotic_burst(x);
      }
    }

    double yL = 0.0;
    double yR = 0.0;

    if (ds.active && ds.burst_total_samples > 0.0) {
      const double env = aureonoise_dichotic_envelope(ds);
      double stim = 0.0;
      if (ds.current_content == t_aureonoise::DichoticContent::Tone) {
        stim = std::sin(ds.tone_phase);
        ds.tone_phase += ds.tone_inc;
        if (ds.tone_phase > aureo::kTwoPi) ds.tone_phase -= aureo::kTwoPi;
      } else {
        ds.noise_state = 0.96 * ds.noise_state + 0.04 * x->rng.uniPM1();
        stim = ds.noise_state;
      }
      const double v = ds.amp * env * stim;
      if (ds.left_ear) yL += v; else yR += v;
      ds.burst_elapsed_samples += 1.0;
      if (ds.burst_elapsed_samples >= ds.burst_total_samples) {
        ds.active = false;
        ds.burst_elapsed_samples = 0.0;
        ds.tone_phase = 0.0;
        ds.noise_state = 0.0;
        ds.samples_until_next = ds.interval_samples;
      }
    }

    outL[n] = yL;
    outR[n] = yR;
    ++x->sample_counter;
  }
}

static bool aureonoise_poisson_enforce(const t_aureonoise* x, double& pan, double minPan, double minSamples)
{
  pan = aureo::clamp(pan, -1.0, 1.0);
  minPan = std::max(0.0, minPan);
  minSamples = std::max(0.0, minSamples);
  if (minPan <= 0.0 || x->grains.empty()) return true;

  double guard = minPan;
  const double relax = 0.82;
  for (int attempt = 0; attempt < 8; ++attempt) {
    bool conflict = false;
    for (const auto& g : x->grains) {
      if (!g.on) continue;
      if (minSamples > 0.0 && static_cast<double>(g.age) > minSamples) continue;
      const double dist = std::abs(pan - g.pan);
      if (dist < guard) {
        conflict = true;
        const double sign = (pan >= g.pan) ? 1.0 : -1.0;
        pan = aureo::clamp(g.pan + sign * guard, -1.0, 1.0);
        break;
      }
    }
    if (!conflict) return true;
    guard *= relax;
  }

  for (const auto& g : x->grains) {
    if (!g.on) continue;
    if (minSamples > 0.0 && static_cast<double>(g.age) > minSamples) continue;
    if (std::abs(pan - g.pan) < 0.5 * minPan) return false;
  }

  pan = aureo::clamp(pan, -1.0, 1.0);
  return true;
}

void* aureonoise_new(t_symbol* s, long argc, t_atom* argv);
void  aureonoise_free(t_aureonoise* x);
void  aureonoise_assist(t_aureonoise* x, void* b, long m, long a, char* s);
void  aureonoise_clear(t_aureonoise* x);
void  aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short* count, double sr, long n, long flags);
void  aureonoise_perform64(t_aureonoise* x, t_object* dsp64, double** ins, long numins,
                           double** outs, long numouts, long sampleframes, long flags, void* userparam);

static t_class* s_aureonoise_class = nullptr;

static uint32_t map_len_samples(t_aureonoise* x, double u)
{
  const double base = aureo::clamp(x->p_baselen_ms, aureo::kMinBaseLengthMs, 2000.0) * 0.001 * x->sr;
  const double kexp = (2.0 * u - 1.0) * aureo::clamp01(x->p_len_phi);
  double L = base * std::pow(aureo::kPhi, kexp);
  L = aureo::clamp(L, aureo::kMinGrainSamples, x->sr * 4.0);
  return static_cast<uint32_t>(L);
}

static int find_free_grain(t_aureonoise* x)
{
  for (int i = 0; i < aureo::kMaxGrains; ++i) {
    if (!x->grains[i].on) return i;
  }
  return -1;
}

static int schedule_gap_samples(t_aureonoise* x)
{
  double rate = aureo::clamp(x->p_rate, 0.0, aureo::kMaxEventRateHz);
#if AUREO_THERMO_LATTICE
  if (x->p_thermo) {
    const double ur = 0.5 + 0.5 * std::tanh(x->ou_rate.y);
    const double rate_phi = aureo::map_phi_range(std::max(0.001, x->p_rate / aureo::kPhi), x->p_rate * aureo::kPhi, ur);
    rate = aureo::clamp(rate_phi, 0.0, aureo::kMaxEventRateHz);
  }
#endif
  if (rate <= 1e-6) {
    const double fallback = std::max(1.0, x->sr * 0.25);
    return static_cast<int>(fallback);
  }

  const double t = static_cast<double>(x->sample_counter) / x->sr;
  double lambda = rate * (1.0 + 0.2 * std::sin(2.0 * aureo::kPi * (t * (1.0 / aureo::kPhi))));
  lambda = std::max(lambda, 1e-3);
#if AUREO_THERMO_LATTICE && AUREO_BURST_HAWKES
  if (x->p_burst) lambda += 0.3 * x->hawkes.lambda;
#endif

  const double U = std::max(1.0e-12, x->rng.uni01());
  double gap_sec = -std::log(U) / lambda;
  const double base_rate = aureo::clamp(x->p_rate, 0.0, aureo::kMaxEventRateHz);
  const double cap_sec = (base_rate > 1e-6) ? std::min(30.0, 4.0 / base_rate) : 0.25;
  gap_sec = std::min(gap_sec, cap_sec);
  return static_cast<int>(std::max(1.0, std::round(gap_sec * x->sr)));
}

static void reset_sequences(t_aureonoise* x)
{
  x->rng.seed(static_cast<uint64_t>(x->p_seed));
  x->w_phi.x = x->rng.uni01();
  x->w_phi.set_step(aureo::kInvPhi);
#if AUREO_WD_PHI_POWERS
  x->w_s2.x = x->rng.uni01();
  x->w_s2.set_step(1.0 / (aureo::kPhi * aureo::kPhi));
  x->w_pl.x = x->rng.uni01();
  x->w_pl.set_step(1.0 / (aureo::kPhi * aureo::kPhi * aureo::kPhi));
#else
  x->w_s2.x = x->rng.uni01();
  x->w_s2.set_step(aureo::kInvSqrt2);
  x->w_pl.x = x->rng.uni01();
  x->w_pl.set_step(aureo::kInvPlastic);
#endif
  x->noise.color = static_cast<aureo::NoiseColor>(x->p_noise_color);
  x->noise.amount = aureo::clamp01(x->p_color_amt);
  x->noise.reset();
  x->tone.clear();
}

extern "C" int C74_EXPORT main(void)
{
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-function-type-mismatch"
#endif
  t_class* c = class_new("aureonoise~",
                         (method)aureonoise_new,
                         (method)aureonoise_free,
                         (long)sizeof(t_aureonoise),
                         0L, A_GIMME, 0);

  class_addmethod(c, (method)aureonoise_assist, "assist", A_CANT, 0);
  class_addmethod(c, (method)aureonoise_clear,  "clear", 0);
  class_addmethod(c, (method)aureonoise_report, "report", 0);
  class_addmethod(c, (method)aureonoise_dsp64,  "dsp64", A_CANT, 0);
  class_dspinit(c);

  aureonoise_setup_attributes(c);

  class_register(CLASS_BOX, c);
  s_aureonoise_class = c;
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
  return 0;
}

#if AUREO_THERMO_LATTICE
void aureonoise_lattice_safe_resize(t_aureonoise* x, int X, int Y, int Z)
{
  if (!x) return;
  X = std::max(2, X);
  Y = std::max(2, Y);
  Z = std::max(1, Z);
  const size_t n = static_cast<size_t>(X) * static_cast<size_t>(Y) * static_cast<size_t>(Z);
  std::vector<double> newx(n, 0.0);
  std::vector<double> newtmp(n, 0.0);
  if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
  x->lat.X = X;
  x->lat.Y = Y;
  x->lat.Z = Z;
  x->lat.x.swap(newx);
  x->lat.tmp.swap(newtmp);
  x->lat_phase = 0.0;
  x->lat_last_v = 0.0;
  if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
}
#endif

void* aureonoise_new(t_symbol*, long argc, t_atom* argv)
{
  auto* x = static_cast<t_aureonoise*>(object_alloc(s_aureonoise_class));
  if (!x) return nullptr;

  dsp_setup(reinterpret_cast<t_pxobject*>(x), 0);
  // In Max gli outlet vengono creati da destra verso sinistra: creiamo prima
  // l'outlet informativo così rimane a destra, quindi il canale destro e infine
  // il sinistro che deve risultare il primo (più a sinistra) sul box.
  x->out_info = outlet_new(reinterpret_cast<t_object*>(x), NULL);
  outlet_new(reinterpret_cast<t_object*>(x), "signal"); // canale destro
  outlet_new(reinterpret_cast<t_object*>(x), "signal"); // canale sinistro

  x->sr = sys_getsr();
  if (x->sr <= 0) x->sr = 44100.0;

  systhread_mutex_new(&x->report_mu, 0);
  aureonoise_reset_reports(x);

  reset_sequences(x);
  make_small_primes(x);
  x->ring.clear();
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  x->lfo_wow_phase = 0.0;
  x->lfo_flut_phase = 0.0;
  for (auto& g : x->grains) g = {};
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;

  x->pinna.width = x->p_width;
  x->pinna.itd_us = x->p_itd_us;
  x->pinna.ild_db = x->p_ild_db;
  x->pinna_enabled = (x->p_pinna_on != 0);
  aureonoise_update_pinna_state(x);
  x->pinna_mix = x->pinna_mix_target;
  aureonoise_update_pinna_filters(x);
  x->pinna_notchL.clear();
  x->pinna_notchR.clear();
  x->field.temperature = 0.45;

#if AUREO_THERMO_LATTICE
  x->ou_pan.tau = 0.60; x->ou_pan.y = 0.0;
  x->ou_itd.tau = 0.40; x->ou_itd.y = 0.0;
  x->ou_amp.tau = 0.80; x->ou_amp.y = 0.0;
  x->ou_rate.tau = 1.20; x->ou_rate.y = 0.0;
  x->lat.init(static_cast<int>(x->p_lat_x), static_cast<int>(x->p_lat_y), static_cast<int>(x->p_lat_z));
  x->lat.eps = x->p_lat_eps;
  x->lat.gamma = x->p_lat_gamma;
  x->lat.sigma = x->p_lat_sigma;
  x->lat_phase = 0.0;
  systhread_mutex_new(&x->lat_mu, 0);
  x->lat_last_v = 0.0;
#if AUREO_BURST_HAWKES
  x->hawkes.base = 4.0;
  x->hawkes.beta = 30.0;
  x->hawkes.lambda = 0.0;
#endif
#endif

  aureonoise_reset_dichotic_state(x);
  attr_args_process(x, static_cast<short>(argc), argv);
  return x;
}

void aureonoise_free(t_aureonoise* x)
{
  dsp_free(reinterpret_cast<t_pxobject*>(x));
  if (x->report_mu) {
    systhread_mutex_free(x->report_mu);
    x->report_mu = nullptr;
  }
  x->out_info = nullptr;
#if AUREO_THERMO_LATTICE
  if (x->lat_mu) {
    systhread_mutex_free(x->lat_mu);
    x->lat_mu = nullptr;
  }
#endif
}

void aureonoise_assist(t_aureonoise*, void*, long m, long a, char* s)
{
  if (m == ASSIST_INLET)
    snprintf_zero(s, 256, "(no inlet)");
  else
    snprintf_zero(s, 256,
                  (a == 0) ? "Out L (aureo glitch noise)" :
                  (a == 1) ? "Out R (aureo glitch noise)" :
                             "Report (event log, diagnostica stocastica)");
}

void aureonoise_report(t_aureonoise* x)
{
  if (!x || !x->out_info) return;

  std::vector<t_aureonoise::GrainReport> snapshot;
  std::vector<t_aureonoise::DichoticReport> dichotic_snapshot;
  if (x->report_mu) systhread_mutex_lock(x->report_mu);
  const size_t count = x->report_count;
  const size_t cap = x->report_log.size();
  if (cap > 0 && count > 0) {
    snapshot.reserve(count);
    size_t idx = (count < cap) ? ((x->report_head + cap - count) % cap) : x->report_head;
    for (size_t i = 0; i < count; ++i) {
      snapshot.push_back(x->report_log[idx]);
      idx = (idx + 1) % cap;
    }
  }
  const size_t dcount = x->dichotic_report_count;
  const size_t dcap = x->dichotic_report_log.size();
  if (dcap > 0 && dcount > 0) {
    dichotic_snapshot.reserve(dcount);
    size_t idx = (dcount < dcap) ? ((x->dichotic_report_head + dcap - dcount) % dcap) : x->dichotic_report_head;
    for (size_t i = 0; i < dcount; ++i) {
      dichotic_snapshot.push_back(x->dichotic_report_log[idx]);
      idx = (idx + 1) % dcap;
    }
  }
  if (x->report_mu) systhread_mutex_unlock(x->report_mu);

  const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
  double sum_gap = 0.0;
  double sum_gap_sq = 0.0;
  for (const auto& r : snapshot) {
    sum_gap += r.gap_samples;
    sum_gap_sq += r.gap_samples * r.gap_samples;
  }
  const double n = static_cast<double>(snapshot.size());
  const double mean_gap_samples = (n > 0.0) ? (sum_gap / n) : 0.0;
  const double mean_gap_sec = mean_gap_samples / sr;
  double var_gap_samples = 0.0;
  if (n > 0.0) {
    const double mean_sq = mean_gap_samples * mean_gap_samples;
    var_gap_samples = std::max(0.0, (sum_gap_sq / n) - mean_sq);
  }
  const double std_gap_sec = std::sqrt(var_gap_samples) / sr;
  const double mean_rate = (mean_gap_sec > 1.0e-9) ? (1.0 / mean_gap_sec) : 0.0;
  const double gap_cv = (mean_gap_sec > 1.0e-9) ? (std_gap_sec / mean_gap_sec) : 0.0;

  size_t dichotic_match_count = 0;
  for (const auto& d : dichotic_snapshot) {
    if (d.match) ++dichotic_match_count;
  }
  const double dichotic_match_ratio = (!dichotic_snapshot.empty())
                                        ? (static_cast<double>(dichotic_match_count) /
                                           static_cast<double>(dichotic_snapshot.size()))
                                        : 0.0;

  static t_symbol* sym_grain = gensym("grain");
  static t_symbol* sym_summary = gensym("summary");
  static t_symbol* sym_id = gensym("id");
  static t_symbol* sym_time = gensym("time_ms");
  static t_symbol* sym_gap = gensym("gap_ms");
  static t_symbol* sym_dur = gensym("dur_ms");
  static t_symbol* sym_amp = gensym("amp");
  static t_symbol* sym_pan = gensym("pan");
  static t_symbol* sym_itd = gensym("itd_samples");
  static t_symbol* sym_ild = gensym("ild_db");
  static t_symbol* sym_binaural_az = gensym("binaural_azimuth_deg");
  static t_symbol* sym_binaural_focus = gensym("binaural_focus");
  static t_symbol* sym_binaural_cf = gensym("binaural_crossfeed");
  static t_symbol* sym_phi = gensym("phi_u");
  static t_symbol* sym_s2 = gensym("s2_u");
  static t_symbol* sym_pl = gensym("pl_u");
  static t_symbol* sym_rng4 = gensym("rng4_u");
  static t_symbol* sym_rng5 = gensym("rng5_u");
  static t_symbol* sym_rng6 = gensym("rng6_u");
  static t_symbol* sym_lat_u = gensym("lat_u");
  static t_symbol* sym_lat_v = gensym("lat_v");
  static t_symbol* sym_lambda = gensym("hawkes_lambda");
  static t_symbol* sym_burst = gensym("burst");
  static t_symbol* sym_thermo = gensym("thermo");
  static t_symbol* sym_lattice = gensym("lattice");
  static t_symbol* sym_dichotic_trial = gensym("dichotic_trial");
  static t_symbol* sym_ear = gensym("ear");
  static t_symbol* sym_match = gensym("match");
  static t_symbol* sym_content = gensym("content");
  static t_symbol* sym_freq = gensym("freq_hz");
  static t_symbol* sym_left = gensym("L");
  static t_symbol* sym_right = gensym("R");
  static t_symbol* sym_tone = gensym("tone");
  static t_symbol* sym_noise = gensym("noise");

  for (const auto& r : snapshot) {
    t_atom av[40];
    long ac = 0;
    atom_setsym(av + ac++, sym_id); atom_setlong(av + ac++, static_cast<t_atom_long>(r.index));
    atom_setsym(av + ac++, sym_time); atom_setfloat(av + ac++, r.timestamp_sec * 1000.0);
    atom_setsym(av + ac++, sym_gap); atom_setfloat(av + ac++, (r.gap_samples / sr) * 1000.0);
    atom_setsym(av + ac++, sym_dur); atom_setfloat(av + ac++, (r.dur_samples / sr) * 1000.0);
    atom_setsym(av + ac++, sym_amp); atom_setfloat(av + ac++, r.amp);
    atom_setsym(av + ac++, sym_pan); atom_setfloat(av + ac++, r.pan);
    atom_setsym(av + ac++, sym_itd); atom_setfloat(av + ac++, r.itd_samples);
    atom_setsym(av + ac++, sym_ild); atom_setfloat(av + ac++, r.ild_db);
    atom_setsym(av + ac++, sym_binaural_az); atom_setfloat(av + ac++, r.binaural_azimuth_deg);
    atom_setsym(av + ac++, sym_binaural_focus); atom_setfloat(av + ac++, r.binaural_focus);
    atom_setsym(av + ac++, sym_binaural_cf); atom_setfloat(av + ac++, r.binaural_crossfeed);
    atom_setsym(av + ac++, sym_phi); atom_setfloat(av + ac++, r.phi_u);
    atom_setsym(av + ac++, sym_s2); atom_setfloat(av + ac++, r.s2_u);
    atom_setsym(av + ac++, sym_pl); atom_setfloat(av + ac++, r.pl_u);
    atom_setsym(av + ac++, sym_rng4); atom_setfloat(av + ac++, r.rng4);
    atom_setsym(av + ac++, sym_rng5); atom_setfloat(av + ac++, r.rng5);
    atom_setsym(av + ac++, sym_rng6); atom_setfloat(av + ac++, r.rng6);
    atom_setsym(av + ac++, sym_lat_u); atom_setfloat(av + ac++, r.lattice_u);
    atom_setsym(av + ac++, sym_lat_v); atom_setfloat(av + ac++, r.lattice_v);
    atom_setsym(av + ac++, sym_lambda); atom_setfloat(av + ac++, r.hawkes_lambda);
    atom_setsym(av + ac++, sym_burst); atom_setlong(av + ac++, r.burst_on ? 1 : 0);
    atom_setsym(av + ac++, sym_thermo); atom_setlong(av + ac++, r.thermo_on ? 1 : 0);
    atom_setsym(av + ac++, sym_lattice); atom_setlong(av + ac++, r.lattice_on ? 1 : 0);
    outlet_anything(x->out_info, sym_grain, ac, av);
  }

  for (const auto& d : dichotic_snapshot) {
    t_atom dv[24];
    long dc = 0;
    atom_setsym(dv + dc++, sym_id); atom_setlong(dv + dc++, static_cast<t_atom_long>(d.index));
    atom_setsym(dv + dc++, sym_time); atom_setfloat(dv + dc++, d.timestamp_sec * 1000.0);
    atom_setsym(dv + dc++, sym_ear); atom_setsym(dv + dc++, d.left_ear ? sym_left : sym_right);
    atom_setsym(dv + dc++, sym_match); atom_setlong(dv + dc++, d.match ? 1 : 0);
    atom_setsym(dv + dc++, sym_content); atom_setsym(dv + dc++, d.is_tone ? sym_tone : sym_noise);
    atom_setsym(dv + dc++, sym_freq); atom_setfloat(dv + dc++, d.freq_hz);
    atom_setsym(dv + dc++, sym_amp); atom_setfloat(dv + dc++, d.amp);
    outlet_anything(x->out_info, sym_dichotic_trial, dc, dv);
  }

  const uint64_t total_events = x->report_total.load(std::memory_order_relaxed);
  const uint64_t dropped_events = x->report_dropped.load(std::memory_order_relaxed);
  const uint64_t total_trials = x->dichotic_report_total.load(std::memory_order_relaxed);
  const uint64_t dropped_trials = x->dichotic_report_dropped.load(std::memory_order_relaxed);

  t_atom summary[24];
  long sc = 0;
  atom_setsym(summary + sc++, gensym("events_total"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(total_events));
  atom_setsym(summary + sc++, gensym("events_buffered"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(snapshot.size()));
  atom_setsym(summary + sc++, gensym("mean_gap_ms"));
  atom_setfloat(summary + sc++, mean_gap_sec * 1000.0);
  atom_setsym(summary + sc++, gensym("mean_rate_hz"));
  atom_setfloat(summary + sc++, mean_rate);
  atom_setsym(summary + sc++, gensym("gap_cv"));
  atom_setfloat(summary + sc++, gap_cv);
  atom_setsym(summary + sc++, gensym("dropped"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(dropped_events));
  atom_setsym(summary + sc++, gensym("dichotic_trials_total"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(total_trials));
  atom_setsym(summary + sc++, gensym("dichotic_trials_buffered"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(dichotic_snapshot.size()));
  atom_setsym(summary + sc++, gensym("dichotic_match_ratio"));
  atom_setfloat(summary + sc++, dichotic_match_ratio);
  atom_setsym(summary + sc++, gensym("dichotic_dropped"));
  atom_setlong(summary + sc++, static_cast<t_atom_long>(dropped_trials));
  outlet_anything(x->out_info, sym_summary, sc, summary);
}

void aureonoise_clear(t_aureonoise* x)
{
  x->ring.clear();
  for (auto& g : x->grains) g = {};
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  x->noise.reset();
  x->pinna_mix = x->pinna_mix_target;
  x->pinna_notchL.clear();
  x->pinna_notchR.clear();
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
  aureonoise_reset_reports(x);
  aureonoise_reset_dichotic_state(x);
#if AUREO_THERMO_LATTICE
  x->lat_phase = 0.0;
  x->lat_last_v = 0.0;
#endif
}

void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short*, double sr, long, long)
{
  x->sr = (sr > 0 ? sr : 44100.0);
  for (auto& g : x->grains) g = {};
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  aureonoise_reset_reports(x);
  x->pinna.width = x->p_width;
  x->pinna.itd_us = x->p_itd_us;
  x->pinna.ild_db = x->p_ild_db;
  const bool was_on = x->pinna_enabled;
  x->pinna_enabled = (x->p_pinna_on != 0);
  if (was_on != x->pinna_enabled) {
    x->pinna_notchL.clear();
    x->pinna_notchR.clear();
  }
  aureonoise_update_pinna_state(x);
  x->pinna_mix = x->pinna_mix_target;
  aureonoise_update_pinna_filters(x);
  aureonoise_reset_dichotic_state(x);
#if AUREO_THERMO_LATTICE
  aureonoise_lattice_safe_resize(x,
                                 static_cast<int>(x->p_lat_x),
                                 static_cast<int>(x->p_lat_y),
                                 static_cast<int>(x->p_lat_z));
  if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
  x->lat.eps = x->p_lat_eps;
  x->lat.gamma = x->p_lat_gamma;
  x->lat.sigma = x->p_lat_sigma;
  if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
#endif
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-function-type-mismatch"
#endif
  object_method(dsp64, gensym("dsp_add64"), x, reinterpret_cast<method>(aureonoise_perform64), 0, NULL);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
}

void aureonoise_perform64(t_aureonoise* x, t_object*, double** ins, long numins,
                          double** outs, long numouts, long sampleframes, long, void*)
{
  if (x->ob.z_disabled) return;
  (void)ins; (void)numins;
  double* outL = (numouts >= 1 && outs[0]) ? outs[0] : nullptr;
  double* outR = (numouts >= 2 && outs[1]) ? outs[1] : nullptr;
  if (!outL || !outR) return;

  if (x->mode == t_aureonoise::StimulusMode::Dichotic) {
    aureonoise_perform64_dichotic(x, outL, outR, sampleframes);
    return;
  }

#ifdef __SSE2__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  if (x->pinna_filters_dirty) {
    aureonoise_update_pinna_filters(x);
  }
  const bool pinna_now = (x->p_pinna_on != 0);
  if (pinna_now != x->pinna_enabled) {
    x->pinna_enabled = pinna_now;
    x->pinna_notchL.clear();
    x->pinna_notchR.clear();
  }
  aureonoise_update_pinna_state(x);
  double pinna_mix = x->pinna_mix;
  const double pinna_target = x->pinna_mix_target;
  const double pinna_depth_amt = aureo::clamp01(x->p_pinna_depth / 24.0);
  double mix_denom = x->sr * 0.02;
  if (mix_denom < 1.0) mix_denom = 1.0;
  const double pinna_mix_step = 1.0 / mix_denom;
  const bool pinna_active = x->pinna_enabled;

  uint32_t wi = x->ring.writeIndex;
  x->noise.color = static_cast<aureo::NoiseColor>(x->p_noise_color);
  x->noise.amount = aureo::clamp01(x->p_color_amt);

  const double wowHz = aureo::map_phi_range(0.1, 1.5, aureo::clamp01(x->p_vhs_wow));
  const double fltHz = aureo::map_phi_range(7.0, 12.0, aureo::clamp01(x->p_vhs_flutter));
  const double incWow = wowHz / x->sr;
  const double incFlt = fltHz / x->sr;
#if AUREO_THERMO_LATTICE
  const double lat_inc = aureo::clamp(x->p_lat_rate, 1.0, 2000.0) / x->sr;
#endif
  const double itd_scale = x->p_itd_us * 1.0e-6 * x->sr;
  const double min_pan_norm = aureo::clamp(x->p_spat_min_deg, 0.0, 180.0) / 90.0;
  const double min_time_samples = aureo::clamp(x->p_spat_min_ms, 0.0, 500.0) * 0.001 * x->sr;
  const double ipd_amt_global = aureo::clamp01(x->p_spat_ipd);
  const double shadow_amt_global = aureo::clamp01(x->p_spat_shadow);
  const auto complex_params = x->field.make_complex_params(x->p_env_complex_mag_var,
                                                           x->p_env_complex_phase_var,
                                                           x->p_env_complex_corr,
                                                           x->p_env_complex_tau_ms,
                                                           x->p_env_complex_bias);

  for (long n = 0; n < sampleframes; ++n) {
    ++x->gap_elapsed;
    x->lfo_wow_phase += incWow; if (x->lfo_wow_phase >= 1.0) x->lfo_wow_phase -= 1.0;
    x->lfo_flut_phase += incFlt; if (x->lfo_flut_phase >= 1.0) x->lfo_flut_phase -= 1.0;
    const double wow = std::sin(2.0 * aureo::kPi * x->lfo_wow_phase);
    const double flt = std::sin(2.0 * aureo::kPi * x->lfo_flut_phase);
    const double vhs_mod = 0.5 * wow + 0.5 * flt;

    pinna_mix += (pinna_target - pinna_mix) * pinna_mix_step;
    pinna_mix = aureo::clamp01(pinna_mix);

#if AUREO_THERMO_LATTICE
    x->lat_phase += lat_inc;
    if (x->lat_phase >= 1.0) {
      const int k = static_cast<int>(std::floor(x->lat_phase));
      const double dt = static_cast<double>(k) / std::max(1.0, x->p_lat_rate);
      bool locked = false;
      if (x->lat_mu) locked = (systhread_mutex_trylock(x->lat_mu) == 0);
      if (!x->lat_mu || locked) {
        for (int i = 0; i < k; ++i) x->lat.step(x->rng);
      }
      if (locked && x->lat_mu) systhread_mutex_unlock(x->lat_mu);
      const double T = aureo::clamp01(x->p_T);
      x->ou_pan.sigma = 0.40 * T;
      x->ou_itd.sigma = 0.35 * T;
      x->ou_amp.sigma = 0.30 * T;
      x->ou_rate.sigma = 0.25 * T;
      x->ou_pan.step(dt, 0.0, x->rng);
      x->ou_itd.step(dt, 0.0, x->rng);
      x->ou_amp.step(dt, 0.0, x->rng);
      x->ou_rate.step(dt, 0.0, x->rng);
#if AUREO_BURST_HAWKES
      (void)x->hawkes.tick(dt, x->rng);
#endif
      x->lat_phase -= static_cast<double>(k);
    }
#endif

    double nz = x->noise.process(x->rng);
#ifdef __aarch64__
    nz += 1.0e-18 * x->rng.uniPM1();
#endif
    nz = aureo::soft_tanh(nz * 1.2);
    x->ring.data[wi] = nz;

    if (--x->samples_to_next <= 0) {
      const int gi = find_free_grain(x);
      bool spawned = false;
      if (gi >= 0) {
        auto& g = x->grains[gi];
        g = {};

        const double u1 = x->w_phi.next();
        const double u2 = x->w_s2.next();
        const double u3 = x->w_pl.next();
        const double u4 = x->rng.uni01();
        const double u5 = x->rng.uni01();
        const double u6 = x->rng.uni01();

        const double gap_samples = static_cast<double>(x->gap_elapsed);
        const double prev_dur = std::max(1.0, x->last_dur_samples);
        const double ratio_gap = gap_samples / prev_dur;
        const double coupling = aureo::clamp01(x->p_hemis_coupling);
        const double time_weight = aureo::clamp(0.35 + 0.65 * (ratio_gap / (ratio_gap + 1.0)), 0.35, 1.0);
        const double hemi = coupling * time_weight;

#if AUREO_THERMO_LATTICE
        double lat_u = 0.5;
        if (x->p_lattice) {
          double v = x->lat_last_v;
          bool locked_lat = false;
          if (x->lat_mu) locked_lat = (systhread_mutex_trylock(x->lat_mu) == 0);
          if (!x->lat_mu || locked_lat) {
            v = x->lat.probe(x->w_phi.next());
            x->lat_last_v = v;
          }
          if (locked_lat && x->lat_mu) systhread_mutex_unlock(x->lat_mu);
          lat_u = 0.5 + 0.5 * std::tanh(v);
        }
        const double oup = x->p_thermo ? aureo::clamp(x->ou_pan.y, -1.0, 1.0) : 0.0;
        const double oua = x->p_thermo ? aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, 0.5 + 0.5 * std::tanh(x->ou_amp.y)) : 1.0;
        const double oui = x->p_thermo ? x->ou_itd.y : 0.0;
#else
        const double lat_u = 0.5;
        const double oup = 0.0;
        const double oua = 1.0;
        const double oui = 0.0;
#endif

        const double amp_shape = std::pow(std::max(1e-9, u2), 0.35);
        const double amp_lat = aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, lat_u);
        const double amp = aureo::kAmpNorm * amp_shape * amp_lat * oua;

        double pan = 2.0 * u3 - 1.0;
#if AUREO_THERMO_LATTICE
        pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
#endif
        pan = aureo::clamp((1.0 - hemi) * pan - hemi * x->prev_pan, -1.0, 1.0);

        bool ok = aureonoise_poisson_enforce(x, pan, min_pan_norm, min_time_samples);
        if (ok) {
          const uint64_t report_idx = x->report_total.fetch_add(1, std::memory_order_relaxed) + 1;
          const double timestamp_sec = (x->sr > 0.0) ? (static_cast<double>(x->sample_counter) / x->sr)
                                                     : 0.0;

          g.dur = map_len_samples(x, u1);
          g.amp = amp;
          g.pan = pan;

          const auto binaural = x->pinna.coefficients(x->sr, pan, u2, u5);
          g.panL = binaural.gainL;
          g.panR = binaural.gainR;
          g.crossfeed = binaural.crossfeed;
          g.focus = binaural.focus;

          double itd = binaural.itd_samples;
#if AUREO_THERMO_LATTICE
          itd += ((2.0 * lat_u - 1.0) * 0.33 + 0.33 * oui) * (x->p_itd_us * 1.0e-6 * x->sr);
#endif
          g.itd = (1.0 - hemi) * itd - hemi * x->prev_itd;

          const double ild_max = aureo::clamp(x->p_ild_db, 0.0, 24.0);
          double ild_db = binaural.ild_db;
          ild_db = aureo::clamp((1.0 - hemi) * ild_db - hemi * x->prev_ild, -ild_max, ild_max);
          g.gL = aureo::db_to_lin(ild_db);
          g.gR = aureo::db_to_lin(-ild_db);

          double ipd_coeff = 0.0;
          if (ipd_amt_global > 1.0e-6) {
            const double ipd_shape = std::abs(2.0 * u6 - 1.0);
            const double ipd_base = 0.18 + 0.55 * ipd_amt_global;
            const double ipd_spread = 0.25 * ipd_amt_global;
            ipd_coeff = aureo::clamp(ipd_base + ipd_spread * (ipd_shape - 0.5), 0.0, 0.95);
          }
          ipd_coeff *= (0.7 + 0.3 * g.focus);
          g.ipd_coeff = ipd_coeff;

          const double shadow_factor = shadow_amt_global * std::abs(pan);
          if (shadow_factor > 1.0e-6) {
            const double fc_min = 800.0;
            const double fc_max = std::max(fc_min, std::min(8000.0, 0.45 * x->sr));
            const double fc = aureo::map_phi_range(fc_min, fc_max, 1.0 - shadow_factor);
            const double alpha = std::exp(-2.0 * aureo::kPi * fc / x->sr);
            g.shadow_a = aureo::clamp(alpha, 0.0, 0.9999);
            g.shadow_left = (pan > 0.0);
            g.shadow_right = (pan < 0.0);
          }

          g.kind = aureo::choose_kind(x->p_glitch_mix, u4);

          int srN = aureo::map_sr_hold_base(x->p_srcrush_amt, u1);
#if AUREO_SR_PRIME_SNAP
          g.sr_holdN = pick_prime_in_range(x, std::max(1, srN - 7), srN + 7, u2);
#else
          g.sr_holdN = srN;
#endif
          g.sr_holdCnt = g.sr_holdN;
          g.q_levels = (1 << (aureo::map_bits(x->p_bitcrush_amt) - 1)) - 1;

          const auto envShape = x->field.make_envelope(x->p_env_attack,
                                                       x->p_env_decay,
                                                       x->p_env_sustain,
                                                       x->p_env_release,
                                                       gap_samples,
                                                       static_cast<double>(g.dur),
                                                       std::abs(pan));
          g.env = envShape;
          g.complex_env = x->field.prepare_complex_envelope(g.dur,
                                                            x->sr,
                                                            x->rng,
                                                            g.env,
                                                            complex_params);
          g.on = true;
          g.age = 0;
          x->prev_pan = pan;
          x->prev_itd = g.itd;
          x->prev_ild = ild_db;
          x->last_gap_samples = gap_samples;
          x->last_dur_samples = static_cast<double>(g.dur);
          x->gap_elapsed = 0;
          t_aureonoise::GrainReport report;
          report.index = report_idx;
          report.timestamp_sec = timestamp_sec;
          report.gap_samples = gap_samples;
          report.dur_samples = static_cast<double>(g.dur);
          report.amp = g.amp;
          report.pan = g.pan;
          report.itd_samples = g.itd;
          report.ild_db = ild_db;
          report.binaural_azimuth_deg = binaural.azimuth_rad * (180.0 / aureo::kPi);
          report.binaural_focus = binaural.focus;
          report.binaural_crossfeed = binaural.crossfeed;
          report.phi_u = u1;
          report.s2_u = u2;
          report.pl_u = u3;
          report.rng4 = u4;
          report.rng5 = u5;
          report.rng6 = u6;
#if AUREO_THERMO_LATTICE
          report.lattice_u = lat_u;
          report.lattice_v = x->lat_last_v;
          report.thermo_on = (x->p_thermo != 0);
          report.lattice_on = (x->p_lattice != 0);
#if AUREO_BURST_HAWKES
          report.hawkes_lambda = x->hawkes.lambda;
          report.burst_on = (x->p_burst != 0) && (x->hawkes.lambda > x->hawkes.base + 1.0);
#else
          report.hawkes_lambda = 0.0;
          report.burst_on = false;
#endif
#else
          report.lattice_u = 0.5;
          report.lattice_v = 0.0;
          report.thermo_on = false;
          report.lattice_on = false;
          report.hawkes_lambda = 0.0;
          report.burst_on = false;
#endif
          aureonoise_log_grain_event(x, report);
          spawned = true;
        }
      }
      x->samples_to_next = schedule_gap_samples(x);
      if (!spawned && gi >= 0) {
        x->grains[gi].on = false;
      }
    }

    double yL = 0.0, yR = 0.0;
    for (auto& g : x->grains) {
      if (!g.on) continue;
      if (g.age >= g.dur) { g.on = false; continue; }

      const double phase = static_cast<double>(g.age) / static_cast<double>(g.dur);
      const double env = x->field.advance_complex_envelope(phase, g.env, g.complex_env, x->rng);

      double itd = g.itd + (vhs_mod * 0.25 * itd_scale);
#if AUREO_THERMO_LATTICE
      if (x->p_lattice) itd += 0.18 * std::tanh(x->lat_last_v) * itd_scale;
#endif
      double sL, sR;
      aureo::ring_read_stereo_itd_frac(x->ring, wi, itd, sL, sR);
      sL *= g.panL * g.gL;
      sR *= g.panR * g.gR;

      if (std::abs(g.crossfeed) > 1.0e-6) {
        const double baseL = sL;
        const double baseR = sR;
        sL = baseL + g.crossfeed * baseR;
        sR = baseR + g.crossfeed * baseL;
      }

      if (--g.sr_holdCnt <= 0) {
        g.heldL = sL;
        g.heldR = sR;
        g.sr_holdCnt = g.sr_holdN;
      }
      sL = g.heldL;
      sR = g.heldR;

      if (g.q_levels > 0) {
        sL = std::round(sL * g.q_levels) / static_cast<double>(g.q_levels);
        sR = std::round(sR * g.q_levels) / static_cast<double>(g.q_levels);
      }

      if (g.kind == aureo::GrainKind::VhsDrop) {
        const double att = 0.5 + 0.5 * (1.0 - std::abs(vhs_mod));
        sL *= att; sR *= att;
      } else if (g.kind == aureo::GrainKind::Stutter) {
        if ((g.age & 7) == 0) { sL *= 0.2; sR *= 0.2; }
      }

      if (g.ipd_coeff > 1.0e-6) {
        sL = aureonoise_allpass(sL, g.ipd_coeff, g.ipd_zL);
        sR = aureonoise_allpass(sR, -g.ipd_coeff, g.ipd_zR);
      }

      if (g.shadow_a > 1.0e-6) {
        if (g.shadow_left) {
          sL = aureonoise_shadow_lp(sL, g.shadow_a, g.shadow_zL);
        }
        if (g.shadow_right) {
          sR = aureonoise_shadow_lp(sR, g.shadow_a, g.shadow_zR);
        }
      }

      yL += g.amp * env * sL;
      yR += g.amp * env * sR;
      ++g.age;
    }

    if (pinna_active) {
      const double dryL = yL;
      const double dryR = yR;
      const double depthAmt = pinna_depth_amt;
      double wetL = x->pinna_notchL.proc(dryL);
      double wetR = x->pinna_notchR.proc(dryR);
      wetL = dryL + depthAmt * (wetL - dryL);
      wetR = dryR + depthAmt * (wetR - dryR);
      const double wet = pinna_mix;
      const double dry_w = 1.0 - wet;
      yL = dry_w * dryL + wet * wetL;
      yR = dry_w * dryR + wet * wetR;
    }

    yL = std::tanh(yL * aureo::kOutDrive) / aureo::kOutDrive;
    yR = std::tanh(yR * aureo::kOutDrive) / aureo::kOutDrive;

    outL[n] = yL;
    outR[n] = yR;

    wi = (wi + 1u) & aureo::kRingMask;
    ++x->sample_counter;
  }

  x->pinna_mix = pinna_mix;
  x->ring.writeIndex = wi;
}

