#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "aureo_config.hpp"
#include "aureo_math.hpp"
#include "aureo_rng.hpp"

namespace aureo {

struct OU {
  double y = 0.0;
  double tau = 0.5;
  double sigma = 0.3;

  double step(double dt, double mu, RNG& rng)
  {
    const double a = std::exp(-dt / std::max(1.0e-6, tau));
    const double s = sigma * std::sqrt(std::max(0.0, 1.0 - a * a));
    y = a * y + (1.0 - a) * mu + s * rng.uniPM1();
    return y;
  }
};

struct Lattice {
  int X = 8, Y = 8, Z = 4;
  double eps = 0.18;
  double gamma = 1.4;
  double sigma = 0.06;
  std::vector<double> x;
  std::vector<double> tmp;

  void init(int x_, int y_, int z_)
  {
    X = x_;
    Y = y_;
    Z = z_;
    x.assign(X * Y * Z, 0.0);
    tmp = x;
  }

  inline int idx(int i, int j, int k) const
  {
    i = (i + X) % X;
    j = (j + Y) % Y;
    k = (k + Z) % Z;
    return (k * Y + j) * X + i;
  }

  inline double act(double v) const { return std::tanh(gamma * v); }

  void step(RNG& rng)
  {
    for (int k = 0; k < Z; ++k) {
      for (int j = 0; j < Y; ++j) {
        for (int i = 0; i < X; ++i) {
          const int p = idx(i, j, k);
          const double s = act(x[idx(i + 1, j, k)]) + act(x[idx(i - 1, j, k)]) +
                           act(x[idx(i, j + 1, k)]) + act(x[idx(i, j - 1, k)]) +
                           act(x[idx(i, j, k + 1)]) + act(x[idx(i, j, k - 1)]);
          tmp[p] = (1.0 - eps) * act(x[p]) + (eps / 6.0) * s + sigma * rng.uniPM1();
        }
      }
    }
    x.swap(tmp);
  }

  double probe(double u) const
  {
    if (x.empty()) return 0.0;
    const size_t n = x.size();
    const size_t p = (size_t)std::floor(clamp01(u) * n) % n;
    return x[p];
  }
};

struct Hawkes {
  double lambda = 0.0;
  double base = 4.0;
  double beta = 30.0;

  bool tick(double dt, RNG& rng)
  {
    lambda = base + (lambda - base) * std::exp(-beta * dt);
    const double p = 1.0 - std::exp(-lambda * dt);
    const bool ev = (rng.uni01() < p);
    if (ev) lambda += base * 0.7;
    return ev;
  }
};

} // namespace aureo

