// aureonoise_tilde.cpp — generatore di texture di rumore/glitch con motore generativo Aureo
// Tutto è spinto da relazioni non-periodiche: φ (rapporto aureo), √2, costante plastica, primi.
// Stereo 0in/2out, pan equal-power + ITD/ILD microbinaurale, “VHS-ish” wow/flutter, dropouts, bit/sr crush.
// Eventi granulari non ripetitivi (sequenze a bassa discrepanza tipo Weyl φ) per ridurre l’assuefazione.
// © 2025 — MIT License. Compilare C++17. Max SDK 8+.
// Nome oggetto in Max: aureonoise~

#include "aureonoise_common.h"

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

// ===================== REGISTRAZIONE CLASSE =====================
static t_class* s_aureonoise_class = nullptr;

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

// ===================== COSTRUZIONE / DISTRUZIONE =====================
void* aureonoise_new(t_symbol* s, long argc, t_atom* argv)
{
  auto* x = (t_aureonoise*)object_alloc(s_aureonoise_class);
  if (!x) return nullptr;

  dsp_setup((t_pxobject*)x, 0);
  outlet_new((t_object*)x, "signal");
  outlet_new((t_object*)x, "signal");

  x->sr = sys_getsr(); if (x->sr <= 0) x->sr = 44100.0;

  // default
  x->p_rate        = 8.0;
  x->p_baselen_ms  = 120.0;
  x->p_len_phi     = 0.8;
  x->p_width       = 1.0;
  x->p_itd_us      = 600.0;
  x->p_ild_db      = 6.0;
  x->p_noise_color = (long)noise_color::pink;
  x->p_color_amt   = 0.65;
  x->p_vhs_wow     = 0.35;
  x->p_vhs_flutter = 0.25;
  x->p_glitch_mix  = 0.5;
  x->p_srcrush_amt = 0.2;
  x->p_bitcrush_amt= 0.15;
  x->p_seed        = 20251010;

  x->rng.seed((uint64_t)x->p_seed);
  make_small_primes(x);

  // Sequenze Weyl
  x->w_phi.x = x->rng.uni01(); x->w_phi.a = kInvPhi;
#if AUREO_WD_PHI_POWERS
  x->w_s2 .x = x->rng.uni01(); x->w_s2 .a = 1.0 / (kPhi*kPhi);
  x->w_pl .x = x->rng.uni01(); x->w_pl .a = 1.0 / (kPhi*kPhi*kPhi);
#else
  x->w_s2 .x = x->rng.uni01(); x->w_s2 .a = kInvSqrt2;
  x->w_pl .x = x->rng.uni01(); x->w_pl .a = kInvPlastic;
#endif

  std::fill(std::begin(x->ring), std::end(x->ring), 0.0);
  x->wi = 0u; x->nz_z1 = x->nz_z2 = 0.0;

  x->lfo_wow_phase = 0.0; x->lfo_flut_phase = 0.0;

  x->sample_counter = 0; x->samples_to_next = (int)(x->sr * 0.1);
  for (auto& g : x->grains) g.on = false;

#if AUREO_THERMO_LATTICE
  // Default Thermo-Lattice
  x->p_thermo    = 1;
  x->p_lattice   = 1;
#if AUREO_BURST_HAWKES
  x->p_burst     = 1;
#else
  x->p_burst     = 0;
#endif
  x->p_T         = 0.45;
  x->p_lat_rate  = 250.0;
  x->p_lat_eps   = 0.18;
  x->p_lat_gamma = 1.40;
  x->p_lat_sigma = 0.06;
  x->p_lat_x = 8; x->p_lat_y = 8; x->p_lat_z = 4;

  x->ou_pan.tau = 0.60; x->ou_itd.tau = 0.40; x->ou_amp.tau = 0.80; x->ou_rate.tau = 1.20;
  x->ou_pan.y = x->ou_itd.y = x->ou_amp.y = x->ou_rate.y = 0.0;
  x->lat.init((int)x->p_lat_x, (int)x->p_lat_y, (int)x->p_lat_z);
  x->lat.eps = x->p_lat_eps; x->lat.gamma = x->p_lat_gamma; x->lat.sigma = x->p_lat_sigma;
  x->lat_phase = 0.0;
#if AUREO_BURST_HAWKES
  x->hawkes.base = 4.0; x->hawkes.beta = 30.0; x->hawkes.lambda = 0.0;
#endif
#endif

  attr_args_process(x, (short)argc, argv);
  return x;
}

void aureonoise_free(t_aureonoise* x) { dsp_free((t_pxobject*)x); }

void aureonoise_assist(t_aureonoise* x, void* b, long m, long a, char* s)
{ if (m == ASSIST_INLET) snprintf_zero(s, 256, "(no inlet)"); else snprintf_zero(s, 256, (a==0) ? "Out L (aureo glitch noise)" : "Out R (aureo glitch noise)"); }

// ===================== COMANDI =====================
void aureonoise_clear(t_aureonoise* x)
{ std::fill(std::begin(x->ring), std::end(x->ring), 0.0); x->wi=0; for(auto& g:x->grains) g.on=false; x->nz_z1=x->nz_z2=0.0; }

// ===================== DSP =====================
void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short*, double sr, long, long)
{
  x->sr = (sr > 0 ? sr : 44100.0);
#if AUREO_THERMO_LATTICE
  x->lat_phase = 0.0;
#endif
  object_method(dsp64, gensym("dsp_add64"), x, (method)aureonoise_perform64, 0, NULL);
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

  uint32_t wi = x->wi;

  // LFO “VHS”: frequenze φ-based (per distribuzioni piacevoli sui range)
  const double wowHz = map_phi_range(0.1, 1.5, clamp01(x->p_vhs_wow));
  const double fltHz = map_phi_range(7.0, 12.0, clamp01(x->p_vhs_flutter));
  const double incWow = wowHz / x->sr;
  const double incFlt = fltHz / x->sr;
#if AUREO_THERMO_LATTICE
  const double lat_inc = clamp(x->p_lat_rate, 1.0, 2000.0) / x->sr;
#endif

  for (long n = 0; n < sampleframes; ++n) {
    // Aggiorna LFO
    x->lfo_wow_phase += incWow; if (x->lfo_wow_phase >= 1.0) x->lfo_wow_phase -= 1.0;
    x->lfo_flut_phase+= incFlt; if (x->lfo_flut_phase>= 1.0) x->lfo_flut_phase-= 1.0;
    double wow = std::sin(2.0*M_PI*x->lfo_wow_phase);
    double flt = std::sin(2.0*M_PI*x->lfo_flut_phase);
    double vhs_mod = 0.5*wow + 0.5*flt; // -1..+1

#if AUREO_THERMO_LATTICE
    // Tick lattice/OU a control-rate
    x->lat_phase += lat_inc;
    if (x->lat_phase >= 1.0) {
      int k = (int)std::floor(x->lat_phase);
      double dt = (double)k / std::max(1.0, x->p_lat_rate);
      for (int i=0;i<k;++i) x->lat.step(x->rng);
      double T = clamp01(x->p_T);
      x->ou_pan.sigma  = 0.40 * T;  x->ou_itd.sigma  = 0.35 * T;
      x->ou_amp.sigma  = 0.30 * T;  x->ou_rate.sigma = 0.25 * T;
      x->ou_pan.step(dt, 0.0, x->rng);
      x->ou_itd.step(dt, 0.0, x->rng);
      x->ou_amp.step(dt, 0.0, x->rng);
      x->ou_rate.step(dt, 0.0, x->rng);
#if AUREO_BURST_HAWKES
      (void)x->hawkes.tick(dt, x->rng);
#endif
      x->lat_phase -= (double)k;
    }
#endif

    // Sorgente di rumore comune
    double nz = colored_noise_sample(x);
    nz = soft_tanh(nz * 1.2); // compressione dolce
    x->ring[wi] = nz;

    // scheduler
    if (--x->samples_to_next <= 0) {
      int gi = find_free_grain(x);
      if (gi >= 0) {
        auto& g = x->grains[gi];
        g.on  = true; g.age = 0;

        // Weyl tridimensionale
        double u1 = x->w_phi.next();
        double u2 = x->w_s2 .next();
        double u3 = x->w_pl .next();
        double u4 = x->rng.uni01();

#if AUREO_THERMO_LATTICE
        double lat_u = 0.5;
        if (x->p_lattice)  lat_u = 0.5 + 0.5 * std::tanh(x->lat.probe(x->w_phi.next()));
        double oup = x->p_thermo ? clamp(x->ou_pan.y, -1.0, 1.0) : 0.0;
        double oua = x->p_thermo ? map_phi_range(1.0/kPhi, kPhi, 0.5 + 0.5*std::tanh(x->ou_amp.y)) : 1.0;
        double oui = x->p_thermo ? x->ou_itd.y : 0.0;
#else
        double lat_u = 0.5; double oup=0.0, oua=1.0, oui=0.0;
#endif

        g.dur = map_len_samples(x, u1);
        g.amp = kAmpNorm * std::pow(std::max(1e-9, u2), 0.35)
                         * map_phi_range(1.0/kPhi, kPhi, lat_u) * oua; // φ @ macro + OU amp

        // pan + equal-power
        double pan = 2.0 * u3 - 1.0;
#if AUREO_THERMO_LATTICE
        pan += 0.25*(2.0*lat_u - 1.0) + 0.35*oup; // micro/macro φ + OU
#endif
        pan = clamp(pan, -1.0, 1.0);
        double gL, gR; pan_equal_power(pan, clamp01(x->p_width), gL, gR);
        g.panL = gL; g.panR = gR;

        // ITD/ILD (ITD frazionario)
        g.itd = map_itd_samples_frac(x, u2);
#if AUREO_THERMO_LATTICE
        g.itd += ((2.0*lat_u - 1.0) * 0.33 + 0.33*oui) * (x->p_itd_us * 1.0e-6 * x->sr);
#endif
        double ildU = u3; g.gL = map_ild_gain(x->p_ild_db, ildU, true); g.gR = map_ild_gain(x->p_ild_db, ildU, false);

        // glitch mode
        g.kind = choose_kind(x, u4);

        int srN = map_sr_hold_base(x->p_srcrush_amt, u1);
#if AUREO_SR_PRIME_SNAP
        g.sr_holdN = pick_prime_in_range(x, std::max(1, srN-7), srN+7, u2);
#else
        g.sr_holdN = srN;
#endif
        g.sr_holdCnt= g.sr_holdN; g.heldL = 0.0; g.heldR = 0.0;

        g.q_levels  = (1 << (map_bits(x->p_bitcrush_amt) - 1)) - 1; // quantizzazione simmetrica
      }
      x->samples_to_next = schedule_gap_samples(x);
    }

    // somma dei grani
    double yL = 0.0, yR = 0.0;
    for (int i=0;i<kMaxGrains;++i) {
      auto& g = x->grains[i]; if (!g.on) continue;
      if (g.age >= g.dur) { g.on = false; continue; }

      double env = hann01((double)g.age / (double)g.dur);

      // ITD micro mod oscillante (VHS) sovrapposto, frazionario
      double itd = g.itd + (vhs_mod * 0.25 * x->p_itd_us * 1.0e-6 * x->sr);

      // lettura coerente dall’anello (ritardo comune)
      double sL, sR;
      ring_read_stereo_itd_frac(x->ring, wi, itd, sL, sR);

      // equal-power pan + ILD
      sL *= g.panL * g.gL;
      sR *= g.panR * g.gR;

      // SR-hold + bitcrush
      if (--g.sr_holdCnt <= 0) {
        g.heldL = sL; g.heldR = sR; g.sr_holdCnt = g.sr_holdN; 
      }
      sL = g.heldL; sR = g.heldR;

      if (g.q_levels > 0) {
        sL = std::round(sL * g.q_levels) / (double)g.q_levels;
        sR = std::round(sR * g.q_levels) / (double)g.q_levels;
      }

      // modalità speciali
      if (g.kind == grain_kind::vhs_drop) {
        double att = 0.5 + 0.5 * (1.0 - std::abs(vhs_mod)); // 0.5..1
        sL *= att; sR *= att;
      } else if (g.kind == grain_kind::stutter) {
        if ((g.age & 7) == 0) { sL *= 0.2; sR *= 0.2; }
      } // aliaser già realizzato con SR/bit crush

      yL += g.amp * env * sL;
      yR += g.amp * env * sR;
      ++g.age;
    }

    // Limiter morbido
    yL = std::tanh(yL * kOutDrive) / kOutDrive;
    yR = std::tanh(yR * kOutDrive) / kOutDrive;

    outL[n] = yL;
    outR[n] = yR;

    wi = (wi + 1u) & kRingMask;
    ++x->sample_counter;
  }
  x->wi = wi;
}
