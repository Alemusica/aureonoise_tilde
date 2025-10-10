#pragma once

#include <cstdint>

namespace aureo {

struct RNG {
  uint64_t s = 0x9E3779B97F4A7C15ULL;
  inline void seed(uint64_t v) { s = v ? v : 0x2545F4914F6CDD1DULL; }
  inline uint64_t nextu64() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * 2685821657736338717ULL;
  }
  inline double uni01()  { return (nextu64() >> 11) * (1.0 / 9007199254740992.0); }
  inline double uniPM1() { return 2.0 * uni01() - 1.0; }
};

} // namespace aureo

