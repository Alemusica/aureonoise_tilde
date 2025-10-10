// aureonoise_tilde.cpp — generatore di texture di rumore/glitch con motore generativo Aureo
// Tutto è spinto da relazioni non‐periodiche: φ (rapporto aureo), √2, costante plastica, primi.
// Stereo 0in/2out, pan equal-power + ITD/ILD microbinaurale, “VHS-ish” wow/flutter, dropouts, bit/sr crush.
// Eventi granulari non ripetitivi (sequenze a bassa discrepanza tipo Weyl φ) per ridurre l’assuefazione.
// © 2025 — MIT License. Compilare C++17. Max SDK 8+.
// Nome oggetto in Max: aureonoise~

// ======= opzioni compile-time =======
#define AUREO_WD_PHI_POWERS  0   // 0: (1/φ, 1/√2, 1/ρ) [default]; 1: (1/φ, 1/φ², 1/φ³)
#define AUREO_ITD_PHI_SHAPE  0   // 0: ITD lineare [-max..+max]; 1: magnitudine ITD φ-shaped heavy-tail
#define AUREO_SR_PRIME_SNAP  1   // 1: SR-hold snappato ai numeri primi (anti-ritmicità)

extern "C" {
  #include "ext.h"
  #include "ext_obex.h"
  #include "z_dsp.h"
}

#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cstring>

#ifdef __SSE2__
  #include <xmmintrin.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ===================== COSTANTI & UTILI =====================
static constexpr double kPhi        = 1.6180339887498948482;           // φ
static constexpr double kInvPhi     = 1.0 / kPhi;                      // 0.618...
static constexpr double kSqrt2      = 1.4142135623730950488;
static constexpr double kInvSqrt2   = 1.0 / kSqrt2;                    // 0.707...
static constexpr double kPlastic    = 1.3247179572447458000;           // ρ
static constexpr double kInvPlastic = 1.0 / kPlastic;                  // 0.755...
static constexpr int    kRingSize   = 131072;                          // power-of-two per ITD/letture frac
static constexpr int    kRingMask   = kRingSize - 1;
static constexpr int    kMaxGrains  = 32;                              // polifonia granulare
static constexpr double kTiny       = 1.0e-30;
static constexpr double kAmpNorm    = 0.55;                            // calibrazione media ampiezza grani
static constexpr double kOutDrive   = 1.2;                             // limiter morbido

static inline double clamp(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
static inline double clamp01(double x) { return clamp(x, 0.0, 1.0); }
static inline double db_to_lin(double dB) { return std::pow(10.0, dB / 20.0); }
static inline double hann01(double p) { p = clamp01(p); return 0.5 - 0.5 * std::cos(2.0 * M_PI * p); }

// RNG: xorshift64* (veloce, deterministico)
struct RNG {
  uint64_t s = 0x9E3779B97F4A7C15ULL; // seed
  inline void seed(uint64_t v) { s = v ? v : 0x2545F4914F6CDD1DULL; }
  inline uint64_t nextu64() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 2685821657736338717ULL; }
  inline double   uni01()  { return (nextu64() >> 11) * (1.0 / 9007199254740992.0); } // [0,1)
  inline double   uniPM1() { return 2.0 * uni01() - 1.0; }                            // [-1,1)
};

// Weyl sequence (bassa discrepanza) su [0,1): x_{n+1} = frac(x_n + α)
struct Weyl {
  double x = 0.123456789;
  double a = kInvPhi; // passo default: 1/φ
  inline void set_step(double step) { a = step; }
  inline double next() { x += a; x -= std::floor(x); return x; }
};

// Colori di rumore semplici (economici ma efficaci)
enum class noise_color : int { white=0, pink, brown };

// Modalità di glitch/evento
enum class grain_kind : int { burst=0, vhs_drop, stutter, aliaser };

// ===================== STRUTTURA OGGETTO =====================
typedef struct _aureonoise {
  t_pxobject  ob;

  // Parametri esposti
  double      p_rate;           // densità eventi [0..50] Hz
  double      p_baselen_ms;     // durata base di un grano [5..2000] ms
  double      p_len_phi;        // ampiezza dell’esponente φ^k su durata [0..2]
  double      p_width;          // larghezza stereo [0..1]
  double      p_itd_us;         // ITD massimo per grano [0..1000] µs
  double      p_ild_db;         // ILD massimo ±dB [0..24]
  long        p_noise_color;    // white/pink/brown
  double      p_color_amt;      // quantità di colorazione [0..1]
  double      p_vhs_wow;        // 0..1 ampiezza wow (0.1–1.5 Hz)
  double      p_vhs_flutter;    // 0..1 ampiezza flutter (7–12 Hz)
  double      p_glitch_mix;     // 0..1 quanto glitchare (probabilità)
  double      p_srcrush_amt;    // 0..1 (aliaser)
  double      p_bitcrush_amt;   // 0..1 (da 16→4 bit)
  long        p_seed;           // seed

  // Stato DSP
  double      sr;
  double      ring[kRingSize];
  uint32_t    wi;               // write index anello rumore comune

  // generativi
  RNG         rng;
  Weyl        w_phi, w_s2, w_pl;

  // Noise color IIR
  double      nz_z1;            // integratore/leaky
  double      nz_z2;            // per brown

  // LFO wow/flutter
  double      lfo_wow_phase;    // cicli
  double      lfo_flut_phase;   // cicli

  // Scheduler
  uint64_t    sample_counter;
  int         samples_to_next;  // countdown verso prossimo avvio grano

  // Primi piccoli
  int         primes[256];
  int         primes_count;

  // Grani attivi
  struct Grain {
    bool on; uint32_t age, dur; double amp;
    double panL, panR; double itd; double gL, gR;
    int sr_holdN, sr_holdCnt; double heldL, heldR;
    int q_levels; grain_kind kind;
  } grains[kMaxGrains];

} t_aureonoise;

// ===================== PROTOTIPI =====================
void make_small_primes(t_aureonoise* x);
int  pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u);
void* aureonoise_new(t_symbol* s, long argc, t_atom* argv);
void  aureonoise_free(t_aureonoise* x);
void  aureonoise_assist(t_aureonoise* x, void* b, long m, long a, char* s);

void  aureonoise_clear(t_aureonoise* x);
void  aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short* count, double sr, long n, long flags);
void  aureonoise_perform64(t_aureonoise* x, t_object* dsp64, double** ins, long numins,
                           double** outs, long numouts, long sampleframes, long flags, void* userparam);

// Attributi setter
t_max_err set_rate(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_baselen(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lenphi(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_width(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_itd(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_ild(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_noise_color(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_color_amt(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_wow(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_flutter(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_glitch_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_srcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_bitcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_seed(t_aureonoise* x, void*, long ac, t_atom* av);

// ===== util auree =====
static inline double map_phi_range(double vmin, double vmax, double u)
{
  u = clamp01(u);
  vmin = std::max(1.0e-12, vmin);
  vmax = std::max(vmin * 1.000001, vmax);
  const double K = std::log(vmax / vmin) / std::log(kPhi);
  return vmin * std::pow(kPhi, K * u);
}

static inline void pan_equal_power(double pan /*-1..+1*/, double width, double& gL, double& gR)
{ double p = clamp(pan * width, -1.0, 1.0); double th = (p + 1.0) * (M_PI * 0.25); gL = std::cos(th); gR = std::sin(th); }

static inline double colored_noise_sample(t_aureonoise* x)
{
  double w = x->rng.uniPM1();
  const noise_color nc = static_cast<noise_color>(x->p_noise_color);
  const double amt = clamp01(x->p_color_amt);
  if (nc == noise_color::white || amt <= 1e-6) return w;

  if (nc == noise_color::pink) {
    x->nz_z1 = (1.0 - 0.02*amt) * x->nz_z1 + (0.02*amt) * w; // ≈ -3 dB/ott
    return (1.0 - amt) * w + amt * x->nz_z1;
  } else { // brown
    x->nz_z2 = (1.0 - 0.001*amt) * x->nz_z2 + (0.03*amt) * w;
    return clamp((1.0 - amt) * w + amt * x->nz_z2, -1.5, 1.5);
  }
}

static inline double soft_tanh(double x) { return std::tanh(x); }

// durata in campioni via φ^k (k∈[-p_len_phi..+p_len_phi])
static inline uint32_t map_len_samples(t_aureonoise* x, double u)
{
  double base = clamp(x->p_baselen_ms, 5.0, 2000.0) * 0.001 * x->sr;
  double kexp = (2.0*u - 1.0) * clamp01(x->p_len_phi);
  double L = base * std::pow(kPhi, kexp);
  L = clamp(L, 32.0, x->sr * 4.0); // 0.7ms .. 4s
  return static_cast<uint32_t>(L);
}

// ITD (campioni): lineare o φ-shaped (heavy-tail). Frazionaria.
static inline double map_itd_samples_frac(t_aureonoise* x, double u)
{
  const double maxSamp = clamp(x->p_itd_us, 0.0, 1000.0) * 1.0e-6 * x->sr;
#if AUREO_ITD_PHI_SHAPE
  // magnitudine da [~0..max], φ-shaped intorno a |2u-1|
  const double m    = std::abs(2.0*u - 1.0);
  const double mag  = map_phi_range(0.0, std::max(0.0, maxSamp), m); // ammetti zero
  const double sign = (u < 0.5) ? -1.0 : +1.0;
  return sign * mag;
#else
  return (2.0*u - 1.0) * maxSamp;
#endif
}

// ILD in lin con magnitudine φ-shaped
static inline double map_ild_gain(double dBmax, double u, bool left)
{
  const double Dmax = clamp(dBmax, 0.0, 24.0);
  const double m    = std::abs(2.0*u - 1.0);
  const double dBmag= map_phi_range(0.0, Dmax, m);
  const double sgn  = (u < 0.5 ? +1.0 : -1.0);
  const double dL   = left ? (sgn * dBmag) : (-sgn * dBmag);
  return db_to_lin(dL);
}

// SR-hold: steps φ-based, con jitter; opzionale snap ai primi
static inline int map_sr_hold_base(double amt, double u)
{
  amt = clamp01(amt);
  double steps = map_phi_range(1.0, 64.0, amt); // φ-based
  int N = (int)std::max(1.0, std::round(steps));
  int jitter = (int)std::floor(clamp01(u) * 0.999 * N);
  return std::max(1, N - jitter);
}

static inline int map_bits(double amt)
{
  amt = clamp01(amt);
  int bits = (int)std::round(16.0 - amt * 12.0);
  return std::max(4, std::min(16, bits));
}

static inline double quantize_bits(double x, int bits)
{ int levels = (1 << (bits - 1)) - 1; return std::round(x * levels) / (double)levels; }

// Scelta tipo di grano
static inline grain_kind choose_kind(t_aureonoise* x, double u)
{
  double m = clamp01(x->p_glitch_mix);
  if (u < 0.25 * m) return grain_kind::vhs_drop;
  if (u < 0.60 * m) return grain_kind::stutter;
  if (u < 1.00 * m) return grain_kind::aliaser;
  return grain_kind::burst;
}

// Prossimo evento (Poisson) con λ modulata
static inline int schedule_gap_samples(t_aureonoise* x)
{
  double rate = clamp(x->p_rate, 0.0, 50.0);
  if (rate <= 1e-6) return (int)(x->sr * 0.25);

  double t = (double)x->sample_counter / x->sr;
  double lambda = rate * (1.0 + 0.2 * std::sin(2.0 * M_PI * (t * (1.0/kPhi))));
  lambda = std::max(lambda, 1e-3);

  double U = std::max(1.0e-12, x->rng.uni01());
  double gap_sec = -std::log(U) / lambda;
  return (int)std::max(1.0, std::round(gap_sec * x->sr));
}

static inline int find_free_grain(t_aureonoise* x)
{ for (int i=0;i<kMaxGrains;++i) if (!x->grains[i].on) return i; return -1; }

// ===== Interpolazione Lagrange 3° (4-tap) =====
static inline double lagrange3(const double* r, double pos)
{
  const int64_t i   = (int64_t)std::floor(pos);
  const double  f   = pos - (double)i;          // 0..1
  const uint32_t i_1= (uint32_t)(i - 1) & kRingMask;
  const uint32_t i0 = (uint32_t)(i + 0) & kRingMask;
  const uint32_t i1 = (uint32_t)(i + 1) & kRingMask;
  const uint32_t i2 = (uint32_t)(i + 2) & kRingMask;

  const double x_1 = r[i_1], x0 = r[i0], x1 = r[i1], x2 = r[i2];
  const double f1  = f - 1.0;
  const double f2  = f - 2.0;

  const double c_1 = -f*(f1)*(f2) / 6.0;
  const double c0  =  (f+1.0)*(f1)*(f2) / 2.0;
  const double c1  = -f*(f+1.0)*(f2) / 2.0;
  const double c2  =  f*(f+1.0)*f1 / 6.0;

  return c_1*x_1 + c0*x0 + c1*x1 + c2*x2;
}

// Lettura stereo con ITD frazionario e ritardo comune (nessun canale anticipa)
static inline void ring_read_stereo_itd_frac(const double* ring, uint32_t base, double d, double& oL, double& oR)
{
  const double ad = std::abs(d);
  const double B  = std::ceil(ad) + 2.0;             // margine per 4 tap
  const double DL = B + (d > 0.0 ? d : 0.0);         // left più in ritardo se d>0
  const double DR = B + (d < 0.0 ? -d : 0.0);        // right più in ritardo se d<0
  const double posL = (double)base - DL;
  const double posR = (double)base - DR;
  oL = lagrange3(ring, posL);
  oR = lagrange3(ring, posR);
}

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

  // ===== Attributi =====
  CLASS_ATTR_DOUBLE(c,  "rate",       0, t_aureonoise, p_rate);
  CLASS_ATTR_ACCESSORS(c, "rate", NULL, set_rate);
  CLASS_ATTR_LABEL(c,  "rate",       0, "Densità eventi (Hz)");
  CLASS_ATTR_SAVE(c,   "rate",       0);

  CLASS_ATTR_DOUBLE(c,  "baselen_ms",0, t_aureonoise, p_baselen_ms);
  CLASS_ATTR_ACCESSORS(c, "baselen_ms", NULL, set_baselen);
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

  CLASS_ATTR_LONG  (c,  "color",     0, t_aureonoise, p_noise_color);
  CLASS_ATTR_ACCESSORS(c, "color", NULL, set_noise_color);
  CLASS_ATTR_STYLE_LABEL(c, "color", 0, "enumindex", "Colore rumore (white/pink/brown)");
#ifdef CLASS_ATTR_ENUMINDEX
  CLASS_ATTR_ENUMINDEX(c, "color", 0, "white pink brown");
#else
  CLASS_ATTR_ENUM(c, "color", 0, "white pink brown");
#endif
  CLASS_ATTR_SAVE(c,   "color",     0);

  CLASS_ATTR_DOUBLE(c,  "color_amt", 0, t_aureonoise, p_color_amt);
  CLASS_ATTR_ACCESSORS(c, "color_amt", NULL, set_color_amt);
  CLASS_ATTR_FILTER_CLIP(c, "color_amt", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "color_amt", 0, "Quantità colorazione");
  CLASS_ATTR_SAVE(c,   "color_amt", 0);

  CLASS_ATTR_DOUBLE(c,  "vhs_wow",   0, t_aureonoise, p_vhs_wow);
  CLASS_ATTR_ACCESSORS(c, "vhs_wow", NULL, set_vhs_wow);
  CLASS_ATTR_FILTER_CLIP(c, "vhs_wow", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "vhs_wow",   0, "Wow (0..1)");
  CLASS_ATTR_SAVE(c,   "vhs_wow",   0);

  CLASS_ATTR_DOUBLE(c,  "vhs_flutter",0, t_aureonoise, p_vhs_flutter);
  CLASS_ATTR_ACCESSORS(c, "vhs_flutter", NULL, set_vhs_flutter);
  CLASS_ATTR_FILTER_CLIP(c, "vhs_flutter", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "vhs_flutter",0, "Flutter (0..1)");
  CLASS_ATTR_SAVE(c,   "vhs_flutter",0);

  CLASS_ATTR_DOUBLE(c,  "glitch_mix",0, t_aureonoise, p_glitch_mix);
  CLASS_ATTR_ACCESSORS(c, "glitch_mix", NULL, set_glitch_mix);
  CLASS_ATTR_FILTER_CLIP(c, "glitch_mix", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "glitch_mix",0, "Probabilità effetti glitch");
  CLASS_ATTR_SAVE(c,   "glitch_mix",0);

  CLASS_ATTR_DOUBLE(c,  "srcrush",   0, t_aureonoise, p_srcrush_amt);
  CLASS_ATTR_ACCESSORS(c, "srcrush", NULL, set_srcrush);
  CLASS_ATTR_FILTER_CLIP(c, "srcrush", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "srcrush",   0, "Sample-Rate crush (0..1)");
  CLASS_ATTR_SAVE(c,   "srcrush",   0);

  CLASS_ATTR_DOUBLE(c,  "bitcrush",  0, t_aureonoise, p_bitcrush_amt);
  CLASS_ATTR_ACCESSORS(c, "bitcrush", NULL, set_bitcrush);
  CLASS_ATTR_FILTER_CLIP(c, "bitcrush", 0.0, 1.0);
  CLASS_ATTR_LABEL(c,  "bitcrush",  0, "Bit crush (0..1)");
  CLASS_ATTR_SAVE(c,   "bitcrush",  0);

  CLASS_ATTR_LONG  (c,  "seed",      0, t_aureonoise, p_seed);
  CLASS_ATTR_ACCESSORS(c, "seed", NULL, set_seed);
  CLASS_ATTR_LABEL(c,  "seed",      0, "Seed");
  CLASS_ATTR_SAVE(c,   "seed",      0);

  class_register(CLASS_BOX, c);
  s_aureonoise_class = c;
  return 0;
}

// ===================== COSTRUZIONE / DISTRUZIONE =====================
void make_small_primes(t_aureonoise* x)
{
  const int limit = 4096; std::vector<bool> comp(limit+1,false); std::vector<int> pl;
  for (int p=2;p*p<=limit;++p) if(!comp[p]) for (int q=p*p;q<=limit;q+=p) comp[q]=true;
  for (int n=2;n<=limit;++n) if(!comp[n]) pl.push_back(n);
  x->primes_count = std::min((int)pl.size(), (int)(sizeof(x->primes)/sizeof(int)));
  for (int i=0;i<x->primes_count;++i) x->primes[i]=pl[i];
}

int pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u)
{
  if (lo>hi) std::swap(lo,hi);
  if (x->primes_count<=0) return std::max(1,lo);
  std::vector<int> cand; cand.reserve(64);
  for (int i=0;i<x->primes_count;++i){ int p=x->primes[i]; if (p>=lo && p<=hi) cand.push_back(p); }
  if (cand.empty()) return std::max(1,lo);
  size_t idx = (size_t)std::floor(clamp01(u) * cand.size());
  if (idx>=cand.size()) idx=cand.size()-1;
  return cand[idx];
}

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

  attr_args_process(x, (short)argc, argv);
  return x;
}

void aureonoise_free(t_aureonoise* x) { dsp_free((t_pxobject*)x); }

void aureonoise_assist(t_aureonoise* x, void* b, long m, long a, char* s)
{ if (m == ASSIST_INLET) snprintf_zero(s, 256, "(no inlet)"); else snprintf_zero(s, 256, (a==0) ? "Out L (aureo glitch noise)" : "Out R (aureo glitch noise)"); }

// ===================== ATTR SETTERS =====================
t_max_err set_rate(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_rate = clamp(atom_getfloat(av),0.0,50.0);} return MAX_ERR_NONE; }
t_max_err set_baselen(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_baselen_ms = clamp(atom_getfloat(av),5.0,2000.0);} return MAX_ERR_NONE; }
t_max_err set_lenphi(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_len_phi = clamp(atom_getfloat(av),0.0,2.0);} return MAX_ERR_NONE; }
t_max_err set_width(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_width = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_itd(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_itd_us = clamp(atom_getfloat(av),0.0,1000.0);} return MAX_ERR_NONE; }
t_max_err set_ild(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_ild_db = clamp(atom_getfloat(av),0.0,24.0);} return MAX_ERR_NONE; }
t_max_err set_noise_color(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ if(atom_gettype(av)==A_SYM){ const char* s=atom_getsym(av)->s_name; if(!std::strcmp(s,"white")) x->p_noise_color=(long)noise_color::white; else if(!std::strcmp(s,"pink")) x->p_noise_color=(long)noise_color::pink; else x->p_noise_color=(long)noise_color::brown; } else { long m=atom_getlong(av); if(m<0)m=0; if(m>2)m=2; x->p_noise_color=m; } } return MAX_ERR_NONE; }
t_max_err set_color_amt(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_color_amt = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_vhs_wow(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_vhs_wow = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_vhs_flutter(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_vhs_flutter = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_glitch_mix(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_glitch_mix = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_srcrush(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_srcrush_amt = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_bitcrush(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_bitcrush_amt = clamp(atom_getfloat(av),0.0,1.0);} return MAX_ERR_NONE; }
t_max_err set_seed(t_aureonoise* x, void*, long ac, t_atom* av){ if(ac&&av){ x->p_seed = (long)atom_getlong(av); x->rng.seed((uint64_t)x->p_seed); } return MAX_ERR_NONE; }

// ===================== COMANDI =====================
void aureonoise_clear(t_aureonoise* x)
{ std::fill(std::begin(x->ring), std::end(x->ring), 0.0); x->wi=0; for(auto& g:x->grains) g.on=false; x->nz_z1=x->nz_z2=0.0; }

// ===================== DSP =====================
void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short*, double sr, long, long)
{
  x->sr = (sr > 0 ? sr : 44100.0);
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

  for (long n = 0; n < sampleframes; ++n) {
    // Aggiorna LFO
    x->lfo_wow_phase += incWow; if (x->lfo_wow_phase >= 1.0) x->lfo_wow_phase -= 1.0;
    x->lfo_flut_phase+= incFlt; if (x->lfo_flut_phase>= 1.0) x->lfo_flut_phase-= 1.0;
    double wow = std::sin(2.0*M_PI*x->lfo_wow_phase);
    double flt = std::sin(2.0*M_PI*x->lfo_flut_phase);
    double vhs_mod = 0.5*wow + 0.5*flt; // -1..+1

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

        g.dur = map_len_samples(x, u1);
        g.amp = kAmpNorm * std::pow(std::max(1e-9, u2), 0.35); // heavy-tail calibrata

        // pan + equal-power
        double pan = 2.0 * u3 - 1.0; double gL, gR; pan_equal_power(pan, clamp01(x->p_width), gL, gR);
        g.panL = gL; g.panR = gR;

        // ITD/ILD (ITD frazionario)
        g.itd = map_itd_samples_frac(x, u2);
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
