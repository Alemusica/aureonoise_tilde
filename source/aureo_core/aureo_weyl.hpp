#pragma once

#include <cmath>

#include "aureo_config.hpp"

namespace aureo {

struct Weyl {
  double x = 0.123456789;
  double a = kInvPhi;
  inline void set_step(double step) { a = step; }
  inline double next() {
    x += a;
    x -= std::floor(x);
    return x;
  }
};

} // namespace aureo

