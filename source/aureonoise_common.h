#pragma once

extern "C" {
  #include "ext.h"
  #include "ext_obex.h"
  #include "z_dsp.h"
}

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef AUREO_WD_PHI_POWERS
#define AUREO_WD_PHI_POWERS 0
#endif

#ifndef AUREO_ITD_PHI_SHAPE
#define AUREO_ITD_PHI_SHAPE 0
#endif

#ifndef AUREO_SR_PRIME_SNAP
#define AUREO_SR_PRIME_SNAP 1
#endif

#ifndef AUREO_THERMO_LATTICE
#define AUREO_THERMO_LATTICE 1
#endif

#ifndef AUREO_BURST_HAWKES
#define AUREO_BURST_HAWKES 1
#endif

static constexpr double kPhi        = 1.6180339887498948482;
static constexpr double kInvPhi     = 1.0 / kPhi;
static constexpr double kSqrt2      = 1.4142135623730950488;
static constexpr double kInvSqrt2   = 1.0 / kSqrt2;
static constexpr double kPlastic    = 1.3247179572447458000;
static constexpr double kInvPlastic = 1.0 / kPlastic;
static constexpr int    kRingSize   = 131072;
static constexpr int    kRingMask   = kRingSize - 1;
static constexpr int    kMaxGrains  = 32;
static constexpr double kTiny       = 1.0e-30;
static constexpr double kAmpNorm    = 0.55;
static constexpr double kOutDrive   = 1.2;

inline double clamp(double x, double a, double b) { return x < a ? a : (x > b ? b : x); }
inline double clamp01(double x) { return clamp(x, 0.0, 1.0); }
inline double db_to_lin(double dB) { return std::pow(10.0, dB / 20.0); }
inline double hann01(double p) { p = clamp01(p); return 0.5 - 0.5 * std::cos(2.0 * M_PI * p); }

struct RNG {
  uint64_t s = 0x9E3779B97F4A7C15ULL;
  inline void seed(uint64_t v) { s = v ? v : 0x2545F4914F6CDD1DULL; }
  inline uint64_t nextu64() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 2685821657736338717ULL; }
  inline double   uni01()  { return (nextu64() >> 11) * (1.0 / 9007199254740992.0); }
  inline double   uniPM1() { return 2.0 * uni01() - 1.0; }
};

struct Weyl {
  double x = 0.123456789;
  double a = kInvPhi;
  inline void set_step(double step) { a = step; }
  inline double next() { x += a; x -= std::floor(x); return x; }
};

enum class noise_color : int { white=0, pink, brown };
enum class grain_kind : int { burst=0, vhs_drop, stutter, aliaser };

#if AUREO_THERMO_LATTICE
struct OU {
  double y=0.0, tau=0.5, sigma=0.3;
  inline double step(double dt, double mu, RNG& rng) {
    double a = std::exp(-dt/tau);
    double s = sigma * std::sqrt(std::max(0.0, 1.0 - a*a));
    y = a*y + (1.0 - a)*mu + s * rng.uniPM1();
    return y;
  }
};

struct Lattice {
  int X=8,Y=8,Z=4; double eps=0.18, gamma=1.4, sigma=0.06;
  std::vector<double> x,tmp;
  void init(int x_,int y_,int z_){ X=x_;Y=y_;Z=z_; x.assign(X*Y*Z,0.0); tmp=x; }
  inline int idx(int i,int j,int k) const { i=(i+X)%X; j=(j+Y)%Y; k=(k+Z)%Z; return (k*Y + j)*X + i; }
  inline double act(double v) const { return std::tanh(gamma*v); }
  void step(RNG& rng){
    for(int k=0;k<Z;++k) for(int j=0;j<Y;++j) for(int i=0;i<X;++i){
      int p=idx(i,j,k);
      double s = act(x[idx(i+1,j,k)])+act(x[idx(i-1,j,k)])+act(x[idx(i,j+1,k)])+act(x[idx(i,j-1,k)])+act(x[idx(i,j,k+1)])+act(x[idx(i,j,k-1)]);
      tmp[p] = (1.0-eps)*act(x[p]) + (eps/6.0)*s + sigma*rng.uniPM1();
    }
    x.swap(tmp);
  }
  double probe(double u) const {
    if (x.empty()) return 0.0;
    size_t n=x.size(); size_t p=(size_t)std::floor(clamp01(u)*n) % n;
    return x[p];
  }
};

#if AUREO_BURST_HAWKES
struct Hawkes {
  double lambda=0.0, base=4.0, beta=30.0;
  inline bool tick(double dt, RNG& rng){
    lambda = base + (lambda - base) * std::exp(-beta*dt);
    double p = 1.0 - std::exp(-lambda*dt);
    bool ev = (rng.uni01() < p);
    if (ev) lambda += base * 0.7;
    return ev;
  }
};
#endif
#endif

struct Grain {
  bool on; uint32_t age, dur; double amp;
  double panL, panR; double itd; double gL, gR;
  int sr_holdN, sr_holdCnt; double heldL, heldR;
  int q_levels; grain_kind kind;
};

typedef struct _aureonoise {
  t_pxobject  ob;
  double      p_rate;
  double      p_baselen_ms;
  double      p_len_phi;
  double      p_width;
  double      p_itd_us;
  double      p_ild_db;
  long        p_noise_color;
  double      p_color_amt;
  double      p_vhs_wow;
  double      p_vhs_flutter;
  double      p_glitch_mix;
  double      p_srcrush_amt;
  double      p_bitcrush_amt;
  long        p_seed;
#if AUREO_THERMO_LATTICE
  long        p_thermo;
  long        p_lattice;
  long        p_burst;
  double      p_T;
  double      p_lat_rate;
  double      p_lat_eps;
  double      p_lat_gamma;
  double      p_lat_sigma;
  long        p_lat_x, p_lat_y, p_lat_z;
#endif
  double      sr;
  double      ring[kRingSize];
  uint32_t    wi;
  RNG         rng;
  Weyl        w_phi, w_s2, w_pl;
  double      nz_z1;
  double      nz_z2;
  double      lfo_wow_phase;
  double      lfo_flut_phase;
  uint64_t    sample_counter;
  int         samples_to_next;
  int         primes[256];
  int         primes_count;
  Grain       grains[kMaxGrains];
#if AUREO_THERMO_LATTICE
  OU          ou_pan, ou_itd, ou_amp, ou_rate;
  Lattice     lat;
  double      lat_phase;
#if AUREO_BURST_HAWKES
  Hawkes      hawkes;
#endif
#endif
} t_aureonoise;

inline double map_phi_range(double vmin, double vmax, double u)
{
  u = clamp01(u);
  vmin = std::max(1.0e-12, vmin);
  vmax = std::max(vmin * 1.000001, vmax);
  const double K = std::log(vmax / vmin) / std::log(kPhi);
  return vmin * std::pow(kPhi, K * u);
}

inline void pan_equal_power(double pan, double width, double& gL, double& gR)
{ double p = clamp(pan * width, -1.0, 1.0); double th = (p + 1.0) * (M_PI * 0.25); gL = std::cos(th); gR = std::sin(th); }

inline double colored_noise_sample(t_aureonoise* x)
{
  double w = x->rng.uniPM1();
  const noise_color nc = static_cast<noise_color>(x->p_noise_color);
  const double amt = clamp01(x->p_color_amt);
  if (nc == noise_color::white || amt <= 1e-6) return w;

  if (nc == noise_color::pink) {
    x->nz_z1 = (1.0 - 0.02*amt) * x->nz_z1 + (0.02*amt) * w;
    return (1.0 - amt) * w + amt * x->nz_z1;
  } else {
    x->nz_z2 = (1.0 - 0.001*amt) * x->nz_z2 + (0.03*amt) * w;
    return clamp((1.0 - amt) * w + amt * x->nz_z2, -1.5, 1.5);
  }
}

inline double soft_tanh(double v) { return std::tanh(v); }

inline uint32_t map_len_samples(t_aureonoise* x, double u)
{
  double base = clamp(x->p_baselen_ms, 5.0, 2000.0) * 0.001 * x->sr;
  double kexp = (2.0*u - 1.0) * clamp01(x->p_len_phi);
  double L = base * std::pow(kPhi, kexp);
  L = clamp(L, 32.0, x->sr * 4.0);
  return static_cast<uint32_t>(L);
}

inline double map_itd_samples_frac(t_aureonoise* x, double u)
{
  const double maxSamp = clamp(x->p_itd_us, 0.0, 1000.0) * 1.0e-6 * x->sr;
#if AUREO_ITD_PHI_SHAPE
  const double m    = std::abs(2.0*u - 1.0);
  const double mag  = map_phi_range(0.0, std::max(0.0, maxSamp), m);
  const double sign = (u < 0.5) ? -1.0 : +1.0;
  return sign * mag;
#else
  return (2.0*u - 1.0) * maxSamp;
#endif
}

inline double map_ild_gain(double dBmax, double u, bool left)
{
  const double Dmax = clamp(dBmax, 0.0, 24.0);
  const double m    = std::abs(2.0*u - 1.0);
  const double dBmag= map_phi_range(0.0, Dmax, m);
  const double sgn  = (u < 0.5 ? +1.0 : -1.0);
  const double dL   = left ? (sgn * dBmag) : (-sgn * dBmag);
  return db_to_lin(dL);
}

inline int map_sr_hold_base(double amt, double u)
{
  amt = clamp01(amt);
  double steps = map_phi_range(1.0, 64.0, amt);
  int N = (int)std::max(1.0, std::round(steps));
  int jitter = (int)std::floor(clamp01(u) * 0.999 * N);
  return std::max(1, N - jitter);
}

inline int map_bits(double amt)
{
  amt = clamp01(amt);
  int bits = (int)std::round(16.0 - amt * 12.0);
  return std::max(4, std::min(16, bits));
}

inline double quantize_bits(double x, int bits)
{ int levels = (1 << (bits - 1)) - 1; return std::round(x * levels) / (double)levels; }

inline grain_kind choose_kind(t_aureonoise* x, double u)
{
  double m = clamp01(x->p_glitch_mix);
  if (u < 0.25 * m) return grain_kind::vhs_drop;
  if (u < 0.60 * m) return grain_kind::stutter;
  if (u < 1.00 * m) return grain_kind::aliaser;
  return grain_kind::burst;
}

inline int schedule_gap_samples(t_aureonoise* x)
{
  double rate = clamp(x->p_rate, 0.0, 50.0);
#if AUREO_THERMO_LATTICE
  if (x->p_thermo) {
    double ur = 0.5 + 0.5 * std::tanh(x->ou_rate.y);
    double rate_phi = map_phi_range(std::max(0.001, x->p_rate / kPhi), x->p_rate * kPhi, ur);
    rate = clamp(rate_phi, 0.0, 50.0);
  }
#endif
  if (rate <= 1e-6) return (int)(x->sr * 0.25);

  double t = (double)x->sample_counter / x->sr;
  double lambda = rate * (1.0 + 0.2 * std::sin(2.0 * M_PI * (t * (1.0/kPhi))));
  lambda = std::max(lambda, 1e-3);
#if AUREO_THERMO_LATTICE && AUREO_BURST_HAWKES
  if (x->p_burst) lambda += 0.3 * x->hawkes.lambda;
#endif

  double U = std::max(1.0e-12, x->rng.uni01());
  double gap_sec = -std::log(U) / lambda;
  return (int)std::max(1.0, std::round(gap_sec * x->sr));
}

inline int find_free_grain(t_aureonoise* x)
{ for (int i=0;i<kMaxGrains;++i) if (!x->grains[i].on) return i; return -1; }

inline double lagrange3(const double* r, double pos)
{
  const int64_t i   = (int64_t)std::floor(pos);
  const double  f   = pos - (double)i;
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

inline void ring_read_stereo_itd_frac(const double* ring, uint32_t base, double d, double& oL, double& oR)
{
  const double ad = std::abs(d);
  const double B  = std::ceil(ad) + 2.0;
  const double DL = B + (d > 0.0 ? d : 0.0);
  const double DR = B + (d < 0.0 ? -d : 0.0);
  const double posL = (double)base - DL;
  const double posR = (double)base - DR;
  oL = lagrange3(ring, posL);
  oR = lagrange3(ring, posR);
}

void make_small_primes(t_aureonoise* x);
int  pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u);
void aureonoise_setup_attributes(t_class* c);

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
