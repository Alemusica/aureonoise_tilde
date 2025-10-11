#include "aureonoise~/aureonoise_state.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

void make_small_primes(t_aureonoise* x)
{
  const int limit = 4096;
  std::vector<bool> comp(limit + 1, false);
  std::vector<int> pl;
  for (int p = 2; p * p <= limit; ++p) {
    if (!comp[p]) {
      for (int q = p * p; q <= limit; q += p) comp[q] = true;
    }
  }
  for (int n = 2; n <= limit; ++n) if (!comp[n]) pl.push_back(n);
  const int maxCount = static_cast<int>(x->primes.size());
  x->primes_count = std::min(static_cast<int>(pl.size()), maxCount);
  for (int i = 0; i < x->primes_count; ++i) x->primes[i] = pl[i];
}

int pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u)
{
  if (lo > hi) std::swap(lo, hi);
  if (x->primes_count <= 0) return std::max(1, lo);
  std::vector<int> cand;
  cand.reserve(64);
  for (int i = 0; i < x->primes_count; ++i) {
    const int p = x->primes[i];
    if (p >= lo && p <= hi) cand.push_back(p);
  }
  if (cand.empty()) return std::max(1, lo);
  size_t idx = static_cast<size_t>(std::floor(aureo::clamp01(u) * cand.size()));
  if (idx >= cand.size()) idx = cand.size() - 1;
  return cand[idx];
}

void aureonoise_setup_attributes(t_class* c)
{
  CLASS_ATTR_DOUBLE(c,  "rate",       0, t_aureonoise, p_rate);
  CLASS_ATTR_ACCESSORS(c, "rate", NULL, set_rate);
  CLASS_ATTR_FILTER_CLIP(c, "rate", 0.0, aureo::kMaxEventRateHz);
  CLASS_ATTR_LABEL(c,  "rate",       0, "Densità eventi (Hz)");
  CLASS_ATTR_SAVE(c,   "rate",       0);

  CLASS_ATTR_DOUBLE(c,  "baselen_ms",0, t_aureonoise, p_baselen_ms);
  CLASS_ATTR_ACCESSORS(c, "baselen_ms", NULL, set_baselen);
  CLASS_ATTR_FILTER_CLIP(c, "baselen_ms", aureo::kMinBaseLengthMs, 2000.0);
  CLASS_ATTR_LABEL(c,  "baselen_ms",0, "Durata base grano (ms)");
  CLASS_ATTR_SAVE(c,   "baselen_ms",0);

  CLASS_ATTR_DOUBLE(c,  "len_phi",   0, t_aureonoise, p_len_phi);
  CLASS_ATTR_ACCESSORS(c, "len_phi", NULL, set_lenphi);
  CLASS_ATTR_FILTER_CLIP(c, "len_phi", 0.0, 2.0);
  CLASS_ATTR_LABEL(c,  "len_phi",   0, "Ampiezza esponente φ su durata");
  CLASS_ATTR_SAVE(c,   "len_phi",   0);

  CLASS_ATTR_DOUBLE(c,  "width",     0, t_aureonoise, p_width);
  CLASS_ATTR_ACCESSORS(c, "width", NULL, set_width);
  CLASS_ATTR_FILTER_CLIP(c, "width", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "width",     0, "Larghezza stereo (0..1)");
  CLASS_ATTR_SAVE(c,   "width",     0);

  CLASS_ATTR_DOUBLE(c,  "itd_us",    0, t_aureonoise, p_itd_us);
  CLASS_ATTR_ACCESSORS(c, "itd_us", NULL, set_itd);
  CLASS_ATTR_LABEL(c,  "itd_us",    0, "ITD max (µs)");
  CLASS_ATTR_SAVE(c,   "itd_us",    0);

  CLASS_ATTR_DOUBLE(c,  "ild_db",    0, t_aureonoise, p_ild_db);
  CLASS_ATTR_ACCESSORS(c, "ild_db", NULL, set_ild);
  CLASS_ATTR_LABEL(c,  "ild_db",    0, "ILD ±dB max");
  CLASS_ATTR_SAVE(c,   "ild_db",    0);

  CLASS_ATTR_DOUBLE(c, "spat_min_deg", 0, t_aureonoise, p_spat_min_deg);
  CLASS_ATTR_ACCESSORS(c, "spat_min_deg", NULL, set_spat_min_deg);
  CLASS_ATTR_FILTER_CLIP(c, "spat_min_deg", 0.0, 180.0);
  CLASS_ATTR_LABEL(c, "spat_min_deg", 0, "Separazione minima angolo (°)");
  CLASS_ATTR_SAVE(c, "spat_min_deg", 0);

  CLASS_ATTR_DOUBLE(c, "spat_min_ms", 0, t_aureonoise, p_spat_min_ms);
  CLASS_ATTR_ACCESSORS(c, "spat_min_ms", NULL, set_spat_min_ms);
  CLASS_ATTR_FILTER_CLIP(c, "spat_min_ms", 0.0, 500.0);
  CLASS_ATTR_LABEL(c, "spat_min_ms", 0, "Separazione minima temporale (ms)");
  CLASS_ATTR_SAVE(c, "spat_min_ms", 0);

  CLASS_ATTR_DOUBLE(c, "spat_ipd", 0, t_aureonoise, p_spat_ipd);
  CLASS_ATTR_ACCESSORS(c, "spat_ipd", NULL, set_spat_ipd);
  CLASS_ATTR_FILTER_CLIP(c, "spat_ipd", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "spat_ipd", 0, "Decorrelazione IPD (0..1)");
  CLASS_ATTR_SAVE(c, "spat_ipd", 0);

  CLASS_ATTR_DOUBLE(c, "spat_shadow", 0, t_aureonoise, p_spat_shadow);
  CLASS_ATTR_ACCESSORS(c, "spat_shadow", NULL, set_spat_shadow);
  CLASS_ATTR_FILTER_CLIP(c, "spat_shadow", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "spat_shadow", 0, "Head-shadow controlaterale (0..1)");
  CLASS_ATTR_SAVE(c, "spat_shadow", 0);

  CLASS_ATTR_DOUBLE(c,  "env_attack", 0, t_aureonoise, p_env_attack);
  CLASS_ATTR_ACCESSORS(c, "env_attack", NULL, set_env_attack);
  CLASS_ATTR_FILTER_CLIP(c, "env_attack", 0.01, 4.0);
  CLASS_ATTR_LABEL(c,  "env_attack", 0, "Rapporto attacco/IOI");
  CLASS_ATTR_SAVE(c,   "env_attack", 0);

  CLASS_ATTR_DOUBLE(c,  "env_decay", 0, t_aureonoise, p_env_decay);
  CLASS_ATTR_ACCESSORS(c, "env_decay", NULL, set_env_decay);
  CLASS_ATTR_FILTER_CLIP(c, "env_decay", 0.01, 4.0);
  CLASS_ATTR_LABEL(c,  "env_decay", 0, "Rapporto decadimento/IOI");
  CLASS_ATTR_SAVE(c,   "env_decay", 0);

  CLASS_ATTR_DOUBLE(c,  "env_sustain", 0, t_aureonoise, p_env_sustain);
  CLASS_ATTR_ACCESSORS(c, "env_sustain", NULL, set_env_sustain);
  CLASS_ATTR_FILTER_CLIP(c, "env_sustain", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "env_sustain", 0, "Livello sustain (0..1)");
  CLASS_ATTR_SAVE(c,   "env_sustain", 0);

  CLASS_ATTR_DOUBLE(c,  "env_release", 0, t_aureonoise, p_env_release);
  CLASS_ATTR_ACCESSORS(c, "env_release", NULL, set_env_release);
  CLASS_ATTR_FILTER_CLIP(c, "env_release", 0.01, 4.0);
  CLASS_ATTR_LABEL(c,  "env_release", 0, "Rapporto release/IOI");
  CLASS_ATTR_SAVE(c,   "env_release", 0);

  CLASS_ATTR_DOUBLE(c,  "hemis_coupling", 0, t_aureonoise, p_hemis_coupling);
  CLASS_ATTR_ACCESSORS(c, "hemis_coupling", NULL, set_hemis_coupling);
  CLASS_ATTR_FILTER_CLIP(c, "hemis_coupling", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "hemis_coupling", 0, "Accoppiamento emisferico (0..1)");
  CLASS_ATTR_SAVE(c,   "hemis_coupling", 0);

  CLASS_ATTR_LONG(c,   "pinna_on",  0, t_aureonoise, p_pinna_on);
  CLASS_ATTR_ACCESSORS(c, "pinna_on", NULL, set_pinna_on);
  CLASS_ATTR_STYLE_LABEL(c, "pinna_on", 0, "onoff", "Pinna notch on/off");
  CLASS_ATTR_SAVE(c,   "pinna_on",  0);

  CLASS_ATTR_DOUBLE(c, "pinna_depth", 0, t_aureonoise, p_pinna_depth);
  CLASS_ATTR_ACCESSORS(c, "pinna_depth", NULL, set_pinna_depth);
  CLASS_ATTR_FILTER_CLIP(c, "pinna_depth", 0.0, 24.0);
  CLASS_ATTR_LABEL(c, "pinna_depth", 0, "Profondità notch pinna (dB)");
  CLASS_ATTR_SAVE(c, "pinna_depth", 0);

  CLASS_ATTR_LONG  (c,  "color",     0, t_aureonoise, p_noise_color);
  CLASS_ATTR_ACCESSORS(c, "color", NULL, set_noise_color);
  CLASS_ATTR_STYLE_LABEL(c, "color", 0, "enumindex", "Colore rumore (white/pink/brown risonante)");
#ifdef CLASS_ATTR_ENUMINDEX
  CLASS_ATTR_ENUMINDEX(c, "color", 0, "white pink brown");
#else
  CLASS_ATTR_ENUM(c, "color", 0, "white pink brown");
#endif
  CLASS_ATTR_SAVE(c,   "color",     0);

  CLASS_ATTR_DOUBLE(c, "color_amt", 0, t_aureonoise, p_color_amt);
  CLASS_ATTR_ACCESSORS(c, "color_amt", NULL, set_color_amt);
  CLASS_ATTR_FILTER_CLIP(c, "color_amt", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "color_amt", 0, "Quantità colorazione");
  CLASS_ATTR_SAVE(c,   "color_amt", 0);

  CLASS_ATTR_DOUBLE(c, "vhs_wow",   0, t_aureonoise, p_vhs_wow);
  CLASS_ATTR_ACCESSORS(c, "vhs_wow", NULL, set_vhs_wow);
  CLASS_ATTR_FILTER_CLIP(c, "vhs_wow", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "vhs_wow",   0, "Wow VHS (0..1)");
  CLASS_ATTR_SAVE(c,   "vhs_wow",   0);

  CLASS_ATTR_DOUBLE(c, "vhs_flutter",0, t_aureonoise, p_vhs_flutter);
  CLASS_ATTR_ACCESSORS(c, "vhs_flutter", NULL, set_vhs_flutter);
  CLASS_ATTR_FILTER_CLIP(c, "vhs_flutter", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "vhs_flutter",0, "Flutter VHS (0..1)");
  CLASS_ATTR_SAVE(c,   "vhs_flutter",0);

  CLASS_ATTR_DOUBLE(c, "glitch_mix",0, t_aureonoise, p_glitch_mix);
  CLASS_ATTR_ACCESSORS(c, "glitch_mix", NULL, set_glitch_mix);
  CLASS_ATTR_FILTER_CLIP(c, "glitch_mix", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "glitch_mix",0, "Mix glitch (0..1)");
  CLASS_ATTR_SAVE(c,   "glitch_mix",0);

  CLASS_ATTR_DOUBLE(c, "srcrush_amt",0, t_aureonoise, p_srcrush_amt);
  CLASS_ATTR_ACCESSORS(c, "srcrush_amt", NULL, set_srcrush);
  CLASS_ATTR_FILTER_CLIP(c, "srcrush_amt", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "srcrush_amt",0, "Sample rate crush");
  CLASS_ATTR_SAVE(c,   "srcrush_amt",0);

  CLASS_ATTR_DOUBLE(c, "bitcrush_amt",0, t_aureonoise, p_bitcrush_amt);
  CLASS_ATTR_ACCESSORS(c, "bitcrush_amt", NULL, set_bitcrush);
  CLASS_ATTR_FILTER_CLIP(c, "bitcrush_amt", 0.0, 1.0);
  CLASS_ATTR_LABEL(c, "bitcrush_amt",0, "Bit crush");
  CLASS_ATTR_SAVE(c,   "bitcrush_amt",0);

  CLASS_ATTR_LONG(c, "seed", 0, t_aureonoise, p_seed);
  CLASS_ATTR_ACCESSORS(c, "seed", NULL, set_seed);
  CLASS_ATTR_LABEL(c, "seed", 0, "Seed RNG");
  CLASS_ATTR_SAVE(c, "seed", 0);

#if AUREO_THERMO_LATTICE
  CLASS_ATTR_LONG(c,   "thermo",    0, t_aureonoise, p_thermo);
  CLASS_ATTR_ACCESSORS(c, "thermo", NULL, set_thermo);
  CLASS_ATTR_STYLE_LABEL(c, "thermo", 0, "onoff", "Thermo OU on/off");
  CLASS_ATTR_SAVE(c,   "thermo",    0);

  CLASS_ATTR_LONG(c,   "lattice",   0, t_aureonoise, p_lattice);
  CLASS_ATTR_ACCESSORS(c, "lattice", NULL, set_lattice);
  CLASS_ATTR_STYLE_LABEL(c, "lattice", 0, "onoff", "Lattice on/off");
  CLASS_ATTR_SAVE(c,   "lattice",   0);

#if AUREO_BURST_HAWKES
  CLASS_ATTR_LONG(c,   "burst",     0, t_aureonoise, p_burst);
  CLASS_ATTR_ACCESSORS(c, "burst", NULL, set_burst);
  CLASS_ATTR_STYLE_LABEL(c, "burst", 0, "onoff", "Burst clustering (Hawkes) on/off");
  CLASS_ATTR_SAVE(c,   "burst",     0);
#endif

  CLASS_ATTR_DOUBLE(c, "T",         0, t_aureonoise, p_T);
  CLASS_ATTR_ACCESSORS(c, "T", NULL, set_T);
  CLASS_ATTR_FILTER_CLIP(c, "T", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "T",         0, "Temperatura stocastica (0..1)");
  CLASS_ATTR_SAVE(c,   "T",         0);

  CLASS_ATTR_DOUBLE(c, "lat_rate",  0, t_aureonoise, p_lat_rate);
  CLASS_ATTR_ACCESSORS(c, "lat_rate", NULL, set_lat_rate);
  CLASS_ATTR_LABEL(c,  "lat_rate",  0, "Frequenza lattice (Hz)");
  CLASS_ATTR_SAVE(c,   "lat_rate",  0);

  CLASS_ATTR_DOUBLE(c, "lat_eps",   0, t_aureonoise, p_lat_eps);
  CLASS_ATTR_ACCESSORS(c, "lat_eps", NULL, set_lat_eps);
  CLASS_ATTR_LABEL(c,  "lat_eps",   0, "Accoppiamento lattice ε");
  CLASS_ATTR_SAVE(c,   "lat_eps",   0);

  CLASS_ATTR_DOUBLE(c, "lat_gamma", 0, t_aureonoise, p_lat_gamma);
  CLASS_ATTR_ACCESSORS(c, "lat_gamma", NULL, set_lat_gamma);
  CLASS_ATTR_LABEL(c,  "lat_gamma", 0, "Nonlinearità γ");
  CLASS_ATTR_SAVE(c,   "lat_gamma", 0);

  CLASS_ATTR_DOUBLE(c, "lat_sigma", 0, t_aureonoise, p_lat_sigma);
  CLASS_ATTR_ACCESSORS(c, "lat_sigma", NULL, set_lat_sigma);
  CLASS_ATTR_LABEL(c,  "lat_sigma", 0, "Rumore lattice σ");
  CLASS_ATTR_SAVE(c,   "lat_sigma", 0);

  CLASS_ATTR_LONG(c,   "lat_x",     0, t_aureonoise, p_lat_x);
  CLASS_ATTR_ACCESSORS(c, "lat_x", NULL, set_lat_x);
  CLASS_ATTR_LABEL(c,  "lat_x",     0, "Lattice X");
  CLASS_ATTR_SAVE(c,   "lat_x",     0);

  CLASS_ATTR_LONG(c,   "lat_y",     0, t_aureonoise, p_lat_y);
  CLASS_ATTR_ACCESSORS(c, "lat_y", NULL, set_lat_y);
  CLASS_ATTR_LABEL(c,  "lat_y",     0, "Lattice Y");
  CLASS_ATTR_SAVE(c,   "lat_y",     0);

  CLASS_ATTR_LONG(c,   "lat_z",     0, t_aureonoise, p_lat_z);
  CLASS_ATTR_ACCESSORS(c, "lat_z", NULL, set_lat_z);
  CLASS_ATTR_LABEL(c,  "lat_z",     0, "Lattice Z");
  CLASS_ATTR_SAVE(c,   "lat_z",     0);
#endif
}

t_max_err set_rate(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_rate = aureo::clamp(atom_getfloat(av), 0.0, aureo::kMaxEventRateHz);
  return MAX_ERR_NONE;
}

t_max_err set_baselen(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_baselen_ms = aureo::clamp(atom_getfloat(av), aureo::kMinBaseLengthMs, 2000.0);
  return MAX_ERR_NONE;
}

t_max_err set_lenphi(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_len_phi = aureo::clamp(atom_getfloat(av), 0.0, 2.0);
  return MAX_ERR_NONE;
}

t_max_err set_width(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_width = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
    x->pinna.width = x->p_width;
    aureonoise_mark_pinna_dirty(x);
    aureonoise_update_pinna_filters(x);
  }
  return MAX_ERR_NONE;
}

t_max_err set_itd(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_itd_us = aureo::clamp(atom_getfloat(av), 0.0, 1000.0);
    x->pinna.itd_us = x->p_itd_us;
  }
  return MAX_ERR_NONE;
}

t_max_err set_ild(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_ild_db = aureo::clamp(atom_getfloat(av), 0.0, 24.0);
    x->pinna.ild_db = x->p_ild_db;
  }
  return MAX_ERR_NONE;
}

t_max_err set_spat_min_deg(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_spat_min_deg = aureo::clamp(atom_getfloat(av), 0.0, 180.0);
  return MAX_ERR_NONE;
}

t_max_err set_spat_min_ms(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_spat_min_ms = aureo::clamp(atom_getfloat(av), 0.0, 500.0);
  return MAX_ERR_NONE;
}

t_max_err set_spat_ipd(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_spat_ipd = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_spat_shadow(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_spat_shadow = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_env_attack(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_env_attack = aureo::clamp(atom_getfloat(av), 0.01, 4.0);
  return MAX_ERR_NONE;
}

t_max_err set_env_decay(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_env_decay = aureo::clamp(atom_getfloat(av), 0.01, 4.0);
  return MAX_ERR_NONE;
}

t_max_err set_env_sustain(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_env_sustain = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_env_release(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_env_release = aureo::clamp(atom_getfloat(av), 0.01, 4.0);
  return MAX_ERR_NONE;
}

t_max_err set_hemis_coupling(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_hemis_coupling = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_noise_color(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    if (atom_gettype(av) == A_SYM) {
      const char* s = atom_getsym(av)->s_name;
      if (!std::strcmp(s, "white"))      x->p_noise_color = static_cast<long>(aureo::NoiseColor::White);
      else if (!std::strcmp(s, "pink")) x->p_noise_color = static_cast<long>(aureo::NoiseColor::Pink);
      else                                x->p_noise_color = static_cast<long>(aureo::NoiseColor::Brown);
    } else {
      long m = atom_getlong(av);
      if (m < 0) m = 0;
      if (m > 2) m = 2;
      x->p_noise_color = m;
    }
    x->noise.color = static_cast<aureo::NoiseColor>(x->p_noise_color);
  }
  return MAX_ERR_NONE;
}

t_max_err set_color_amt(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_color_amt = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
    x->noise.amount = x->p_color_amt;
  }
  return MAX_ERR_NONE;
}

t_max_err set_vhs_wow(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_vhs_wow = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_vhs_flutter(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_vhs_flutter = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_glitch_mix(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_glitch_mix = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_srcrush(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_srcrush_amt = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_bitcrush(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_bitcrush_amt = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
  return MAX_ERR_NONE;
}

t_max_err set_seed(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_seed = static_cast<long>(atom_getlong(av));
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
    x->noise.reset();
  }
  return MAX_ERR_NONE;
}

t_max_err set_pinna_on(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    const long on = atom_getlong(av) ? 1 : 0;
    const bool was_on = (x->p_pinna_on != 0);
    x->p_pinna_on = on;
    aureonoise_update_pinna_state(x);
    const bool now_on = (x->p_pinna_on != 0);
    if (was_on != now_on) {
      x->pinna_notchL.clear();
      x->pinna_notchR.clear();
    }
    x->pinna_enabled = now_on;
    aureonoise_update_pinna_filters(x);
  }
  return MAX_ERR_NONE;
}

t_max_err set_pinna_depth(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_pinna_depth = aureo::clamp(atom_getfloat(av), 0.0, 24.0);
    aureonoise_update_pinna_state(x);
  }
  return MAX_ERR_NONE;
}

#if AUREO_THERMO_LATTICE
t_max_err set_thermo(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_thermo = atom_getlong(av) ? 1 : 0;
  // kick di sicurezza per evitare attese troppo lunghe dopo un toggle
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.02));
  return MAX_ERR_NONE;
}

t_max_err set_lattice(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_lattice = atom_getlong(av) ? 1 : 0;
  return MAX_ERR_NONE;
}

t_max_err set_burst(t_aureonoise* x, void*, long ac, t_atom* av)
{
#if AUREO_BURST_HAWKES
  if (ac && av) x->p_burst = atom_getlong(av) ? 1 : 0;
#else
  (void)x; (void)ac; (void)av;
#endif
  return MAX_ERR_NONE;
}

t_max_err set_T(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_T = aureo::clamp(atom_getfloat(av), 0.0, 1.0);
    x->field.temperature = x->p_T;
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_rate(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) x->p_lat_rate = aureo::clamp(atom_getfloat(av), 1.0, 2000.0);
  return MAX_ERR_NONE;
}

t_max_err set_lat_eps(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_eps = atom_getfloat(av);
    if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
    x->lat.eps = x->p_lat_eps;
    if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_gamma(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_gamma = atom_getfloat(av);
    if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
    x->lat.gamma = x->p_lat_gamma;
    if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_sigma(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_sigma = atom_getfloat(av);
    if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
    x->lat.sigma = x->p_lat_sigma;
    if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_x(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_x = std::max<long>(2, atom_getlong(av));
    aureonoise_lattice_safe_resize(x,
                                   static_cast<int>(x->p_lat_x),
                                   static_cast<int>(x->p_lat_y),
                                   static_cast<int>(x->p_lat_z));
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_y(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_y = std::max<long>(2, atom_getlong(av));
    aureonoise_lattice_safe_resize(x,
                                   static_cast<int>(x->p_lat_x),
                                   static_cast<int>(x->p_lat_y),
                                   static_cast<int>(x->p_lat_z));
  }
  return MAX_ERR_NONE;
}

t_max_err set_lat_z(t_aureonoise* x, void*, long ac, t_atom* av)
{
  if (ac && av) {
    x->p_lat_z = std::max<long>(1, atom_getlong(av));
    aureonoise_lattice_safe_resize(x,
                                   static_cast<int>(x->p_lat_x),
                                   static_cast<int>(x->p_lat_y),
                                   static_cast<int>(x->p_lat_z));
  }
  return MAX_ERR_NONE;
}
#endif

