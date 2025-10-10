#pragma once

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

#include <cstddef>

namespace aureo {
inline constexpr double kPhi        = 1.6180339887498948482;
inline constexpr double kInvPhi     = 1.0 / kPhi;
inline constexpr double kSqrt2      = 1.4142135623730950488;
inline constexpr double kInvSqrt2   = 1.0 / kSqrt2;
inline constexpr double kPlastic    = 1.3247179572447458000;
inline constexpr double kInvPlastic = 1.0 / kPlastic;
inline constexpr double kPi         = 3.14159265358979323846;
inline constexpr int    kRingSize   = 131072;
inline constexpr int    kRingMask   = kRingSize - 1;
inline constexpr int    kMaxGrains  = 32;
inline constexpr double kTiny       = 1.0e-30;
inline constexpr double kAmpNorm    = 0.55;
inline constexpr double kOutDrive   = 1.2;
inline constexpr double kMaxEventRateHz = 240.0;
inline constexpr double kMinBaseLengthMs = 1.0;
inline constexpr double kMinGrainSamples = 8.0;
} // namespace aureo

