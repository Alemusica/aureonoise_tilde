#pragma once

#include "aureo_math.hpp"

namespace aureo {

struct FieldState {
  double temperature = 0.45;

  double envelope(double phase) const
  {
    return hann01(phase);
  }
};

} // namespace aureo

