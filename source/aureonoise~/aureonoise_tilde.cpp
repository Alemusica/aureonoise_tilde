// aureonoise_tilde.cpp — generatore di texture di rumore/glitch con motore generativo Aureo
// Tutto è spinto da relazioni non-periodiche: φ, √2, costante plastica, primi.
// Stereo 0in/2out, pan equal-power + ITD/ILD microbinaurale, “VHS-ish” wow/flutter, dropouts, bit/sr crush.
// Eventi granulari non ripetitivi (sequenze a bassa discrepanza tipo Weyl φ) per ridurre l’assuefazione.
// © 2025 — MIT License. Compilare C++17. Max SDK 8+. Nome oggetto in Max: aureonoise~

#include "aureonoise_state.hpp"

#include <cmath>

#ifdef __SSE2__
#include <xmmintrin.h>
#endif

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
  const double base = aureo::clamp(x->p_baselen_ms, 5.0, 2000.0) * 0.001 * x->sr;
  const double kexp = (2.0 * u - 1.0) * aureo::clamp01(x->p_len_phi);
  double L = base * std::pow(aureo::kPhi, kexp);
  L = aureo::clamp(L, 32.0, x->sr * 4.0);
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
  double rate = aureo::clamp(x->p_rate, 0.0, 50.0);
#if AUREO_THERMO_LATTICE
  if (x->p_thermo) {
    const double ur = 0.5 + 0.5 * std::tanh(x->ou_rate.y);
    const double rate_phi = aureo::map_phi_range(std::max(0.001, x->p_rate / aureo::kPhi), x->p_rate * aureo::kPhi, ur);
    rate = aureo::clamp(rate_phi, 0.0, 50.0);
  }
#endif
  if (rate <= 1e-6) return static_cast<int>(x->sr * 0.25);

  const double t = static_cast<double>(x->sample_counter) / x->sr;
  double lambda = rate * (1.0 + 0.2 * std::sin(2.0 * aureo::kPi * (t * (1.0 / aureo::kPhi))));
  lambda = std::max(lambda, 1e-3);
#if AUREO_THERMO_LATTICE && AUREO_BURST_HAWKES
  if (x->p_burst) lambda += 0.3 * x->hawkes.lambda;
#endif

  const double U = std::max(1.0e-12, x->rng.uni01());
  const double gap_sec = -std::log(U) / lambda;
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
  t_class* c = class_new("aureonoise~",
                         (method)aureonoise_new,
                         (method)aureonoise_free,
                         (long)sizeof(t_aureonoise),
                         0L, A_GIMME, 0);

  class_addmethod(c, (method)aureonoise_assist, "assist", A_CANT, 0);
  class_addmethod(c, (method)aureonoise_clear,  "clear", 0);
  class_addmethod(c, (method)aureonoise_dsp64,  "dsp64", A_CANT, 0);
  class_dspinit(c);

  aureonoise_setup_attributes(c);

  class_register(CLASS_BOX, c);
  s_aureonoise_class = c;
  return 0;
}

void* aureonoise_new(t_symbol*, long argc, t_atom* argv)
{
  auto* x = static_cast<t_aureonoise*>(object_alloc(s_aureonoise_class));
  if (!x) return nullptr;

  dsp_setup(reinterpret_cast<t_pxobject*>(x), 0);
  outlet_new(reinterpret_cast<t_object*>(x), "signal");
  outlet_new(reinterpret_cast<t_object*>(x), "signal");

  x->sr = sys_getsr();
  if (x->sr <= 0) x->sr = 44100.0;

  reset_sequences(x);
  make_small_primes(x);
  x->ring.clear();
  x->samples_to_next = static_cast<int>(x->sr * 0.1);
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
#if AUREO_BURST_HAWKES
  x->hawkes.base = 4.0;
  x->hawkes.beta = 30.0;
  x->hawkes.lambda = 0.0;
#endif
#endif

  attr_args_process(x, static_cast<short>(argc), argv);
  return x;
}

void aureonoise_free(t_aureonoise* x)
{
  dsp_free(reinterpret_cast<t_pxobject*>(x));
}

void aureonoise_assist(t_aureonoise*, void*, long m, long a, char* s)
{
  if (m == ASSIST_INLET)
    snprintf_zero(s, 256, "(no inlet)");
  else
    snprintf_zero(s, 256, (a == 0) ? "Out L (aureo glitch noise)" : "Out R (aureo glitch noise)");
}

void aureonoise_clear(t_aureonoise* x)
{
  x->ring.clear();
  for (auto& g : x->grains) g = {};
  x->samples_to_next = static_cast<int>(x->sr * 0.1);
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
}

void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short*, double sr, long, long)
{
  x->sr = (sr > 0 ? sr : 44100.0);
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
#if AUREO_THERMO_LATTICE
  x->lat_phase = 0.0;
#endif
  x->gap_elapsed = 0;
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
  object_method(dsp64, gensym("dsp_add64"), x, reinterpret_cast<method>(aureonoise_perform64), 0, NULL);
}

void aureonoise_perform64(t_aureonoise* x, t_object*, double** ins, long numins,
                          double** outs, long numouts, long sampleframes, long, void*)
{
  if (x->ob.z_disabled) return;
  (void)ins; (void)numins;
  double* outL = (numouts >= 1 && outs[0]) ? outs[0] : nullptr;
  double* outR = (numouts >= 2 && outs[1]) ? outs[1] : nullptr;
  if (!outL || !outR) return;

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
  const double pinna_depth_norm = aureo::clamp01(x->pinna_depth / 24.0);
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
      for (int i = 0; i < k; ++i) x->lat.step(x->rng);
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
    nz += 1.0e-20;
#endif
    nz = aureo::soft_tanh(nz * 1.2);
    x->ring.data[wi] = nz;

    if (--x->samples_to_next <= 0) {
      const int gi = find_free_grain(x);
      if (gi >= 0) {
        auto& g = x->grains[gi];
        g = {};
        g.on = true;
        g.age = 0;

        const double u1 = x->w_phi.next();
        const double u2 = x->w_s2.next();
        const double u3 = x->w_pl.next();
        const double u4 = x->rng.uni01();
        const double u5 = x->rng.uni01();

        const double gap_samples = static_cast<double>(x->gap_elapsed);
        const double prev_dur = std::max(1.0, x->last_dur_samples);
        const double ratio_gap = gap_samples / prev_dur;
        const double coupling = aureo::clamp01(x->p_hemis_coupling);
        const double time_weight = aureo::clamp(0.35 + 0.65 * (ratio_gap / (ratio_gap + 1.0)), 0.35, 1.0);
        const double hemi = coupling * time_weight;

#if AUREO_THERMO_LATTICE
        double lat_u = 0.5;
        if (x->p_lattice) lat_u = 0.5 + 0.5 * std::tanh(x->lat.probe(x->w_phi.next()));
        const double oup = x->p_thermo ? aureo::clamp(x->ou_pan.y, -1.0, 1.0) : 0.0;
        const double oua = x->p_thermo ? aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, 0.5 + 0.5 * std::tanh(x->ou_amp.y)) : 1.0;
        const double oui = x->p_thermo ? x->ou_itd.y : 0.0;
#else
        const double lat_u = 0.5;
        const double oup = 0.0;
        const double oua = 1.0;
        const double oui = 0.0;
#endif

        g.dur = map_len_samples(x, u1);
        g.amp = aureo::kAmpNorm * std::pow(std::max(1e-9, u2), 0.35) * aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, lat_u) * oua;

        double pan = 2.0 * u3 - 1.0;
#if AUREO_THERMO_LATTICE
        pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
#endif
        pan = aureo::clamp((1.0 - hemi) * pan - hemi * x->prev_pan, -1.0, 1.0);
        auto gains = x->pinna.gains(pan);
        g.panL = gains.first;
        g.panR = gains.second;

        double itd = aureo::map_itd_samples_frac(x->sr, x->p_itd_us, u2);
#if AUREO_THERMO_LATTICE
        itd += ((2.0 * lat_u - 1.0) * 0.33 + 0.33 * oui) * (x->p_itd_us * 1.0e-6 * x->sr);
#endif
        g.itd = (1.0 - hemi) * itd - hemi * x->prev_itd;

        const double ild_max = aureo::clamp(x->p_ild_db, 0.0, 24.0);
        const double ild_sign = 2.0 * u5 - 1.0;
        const double ild_mag = aureo::map_phi_range(0.0, ild_max, std::abs(ild_sign));
        double ild_db = ild_sign * ild_mag;
        ild_db = aureo::clamp((1.0 - hemi) * ild_db - hemi * x->prev_ild, -ild_max, ild_max);
        g.gL = aureo::db_to_lin(ild_db);
        g.gR = aureo::db_to_lin(-ild_db);

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
        x->prev_pan = pan;
        x->prev_itd = g.itd;
        x->prev_ild = ild_db;
        x->last_gap_samples = gap_samples;
        x->last_dur_samples = static_cast<double>(g.dur);
        x->gap_elapsed = 0;
      }
      x->samples_to_next = schedule_gap_samples(x);
    }

    double yL = 0.0, yR = 0.0;
    for (auto& g : x->grains) {
      if (!g.on) continue;
      if (g.age >= g.dur) { g.on = false; continue; }

      const double phase = static_cast<double>(g.age) / static_cast<double>(g.dur);
      const double env = x->field.envelope(phase, g.env);

      double itd = g.itd + (vhs_mod * 0.25 * x->p_itd_us * 1.0e-6 * x->sr);
      double sL, sR;
      aureo::ring_read_stereo_itd_frac(x->ring, wi, itd, sL, sR);
      sL *= g.panL * g.gL;
      sR *= g.panR * g.gR;

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

      yL += g.amp * env * sL;
      yR += g.amp * env * sR;
      ++g.age;
    }

    if (pinna_active) {
      const double w = pinna_mix * pinna_depth_norm;
      const double dryL = yL;
      const double dryR = yR;
      const double wetL = x->pinna_notchL.proc(dryL);
      const double wetR = x->pinna_notchR.proc(dryR);
      const double dry_w = 1.0 - w;
      yL = dry_w * dryL + w * wetL;
      yR = dry_w * dryR + w * wetR;
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

