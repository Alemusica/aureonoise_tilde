#pragma once

extern "C" {
#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
}

#include <array>
#include <cstdint>

#include "aureo_core/aureo_binaural.hpp"
#include "aureo_core/aureo_config.hpp"
#include "aureo_core/aureo_field_adsr.hpp"
#include "aureo_core/aureo_math.hpp"
#include "aureo_core/aureo_noise.hpp"
#include "aureo_core/aureo_pinna.hpp"
#include "aureo_core/aureo_ring.hpp"
#include "aureo_core/aureo_rng.hpp"
#include "aureo_core/aureo_stoch.hpp"
#include "aureo_core/aureo_weyl.hpp"

struct t_aureonoise {
  t_pxobject ob;

  double p_rate = 8.0;
  double p_baselen_ms = 120.0;
  double p_len_phi = 0.8;
  double p_width = 1.0;
  double p_itd_us = 600.0;
  double p_ild_db = 6.0;
  double p_env_attack = 0.18;
  double p_env_decay = 0.28;
  double p_env_sustain = 0.55;
  double p_env_release = 0.30;
  double p_hemis_coupling = 0.6;
  long   p_noise_color = static_cast<long>(aureo::NoiseColor::Pink);
  double p_color_amt = 0.65;
  double p_vhs_wow = 0.35;
  double p_vhs_flutter = 0.25;
  double p_glitch_mix = 0.5;
  double p_srcrush_amt = 0.2;
  double p_bitcrush_amt = 0.15;
  long   p_seed = 20251010;
  long   p_pinna_on = 0;
  double p_pinna_depth = 12.0;

#if AUREO_THERMO_LATTICE
  long p_thermo = 1;
  long p_lattice = 1;
  long p_burst = AUREO_BURST_HAWKES ? 1 : 0;
  double p_T = 0.45;
  double p_lat_rate = 250.0;
  double p_lat_eps = 0.18;
  double p_lat_gamma = 1.40;
  double p_lat_sigma = 0.06;
  long   p_lat_x = 8;
  long   p_lat_y = 8;
  long   p_lat_z = 4;
#endif

  double sr = 44100.0;
  aureo::RNG rng;
  aureo::Weyl w_phi;
  aureo::Weyl w_s2;
  aureo::Weyl w_pl;
  aureo::NoiseColorState noise;
  aureo::FieldState field;
  aureo::RingState ring;
  aureo::Pinna pinna;
  aureo::Biquad tone;
  aureo::Biquad pinna_notchL;
  aureo::Biquad pinna_notchR;

  double pinna_mix = 0.0;
  double pinna_mix_target = 0.0;
  double pinna_depth = 12.0;
  double pinna_freq_left = 7200.0;
  double pinna_freq_right = 8200.0;
  double pinna_q = 5.0;
  double pinna_min_freq = 400.0;
  double pinna_notch_base = 7800.0;
  double pinna_notch_spread = 2200.0;
  bool   pinna_enabled = false;
  bool   pinna_filters_dirty = true;

  double lfo_wow_phase = 0.0;
  double lfo_flut_phase = 0.0;

  uint64_t sample_counter = 0;
  int samples_to_next = 0;
  int gap_elapsed = 0;
  std::array<int, 256> primes{};
  int primes_count = 0;
  std::array<aureo::Grain, aureo::kMaxGrains> grains{};

  double prev_pan = 0.0;
  double prev_itd = 0.0;
  double prev_ild = 0.0;
  double last_gap_samples = 0.0;
  double last_dur_samples = 0.0;

#if AUREO_THERMO_LATTICE
  aureo::OU ou_pan;
  aureo::OU ou_itd;
  aureo::OU ou_amp;
  aureo::OU ou_rate;
  aureo::Lattice lat;
  double lat_phase = 0.0;
#if AUREO_BURST_HAWKES
  aureo::Hawkes hawkes;
#endif
#endif
};

inline void aureonoise_update_pinna_state(t_aureonoise* x)
{
  x->pinna_depth = aureo::clamp(x->p_pinna_depth, 0.0, 24.0);
  x->pinna_mix_target = (x->p_pinna_on != 0) ? 1.0 : 0.0;
}

inline double aureonoise_pinna_weight(const t_aureonoise* x)
{
  return aureo::clamp01(x->pinna_mix) * aureo::clamp01(x->pinna_depth / 24.0);
}

inline void aureonoise_update_pinna_filters(t_aureonoise* x)
{
  double sr = x->sr;
  if (sr <= 0.0) sr = 44100.0;
  if (sr < 4000.0) sr = 4000.0;

  const double width = aureo::clamp(x->p_width, 0.0, 1.0);
  const double base = x->pinna_notch_base;
  const double spread = x->pinna_notch_spread;
  x->pinna_freq_left = aureo::clamp(base - 0.5 * spread * width, 1200.0, 16000.0);
  x->pinna_freq_right = aureo::clamp(base + 0.5 * spread * width, 1200.0, 16000.0);

  x->pinna_notchL.setNotch(sr, x->pinna_freq_left, x->pinna_q, x->pinna_min_freq);
  x->pinna_notchR.setNotch(sr, x->pinna_freq_right, x->pinna_q, x->pinna_min_freq);
  x->pinna_filters_dirty = false;
}

inline void aureonoise_mark_pinna_dirty(t_aureonoise* x)
{
  x->pinna_filters_dirty = true;
}

void aureonoise_setup_attributes(t_class* c);
int  pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u);
void make_small_primes(t_aureonoise* x);

t_max_err set_rate(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_baselen(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lenphi(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_width(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_itd(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_ild(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_attack(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_decay(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_sustain(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_release(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_hemis_coupling(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_noise_color(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_color_amt(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_wow(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_flutter(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_glitch_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_srcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_bitcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_seed(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_pinna_on(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_pinna_depth(t_aureonoise* x, void*, long ac, t_atom* av);
#if AUREO_THERMO_LATTICE
t_max_err set_thermo(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lattice(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_burst(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_T(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_rate(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_eps(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_gamma(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_sigma(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_x(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_y(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_z(t_aureonoise* x, void*, long ac, t_atom* av);
#endif

