//! aureonoise - DSP core constants
//! All constants based on φ (golden ratio) and related mathematical constants

/// φ (phi) - Golden ratio
pub const PHI: f64 = 1.6180339887498948482;

/// 1/φ - Inverse golden ratio
pub const INV_PHI: f64 = 0.6180339887498948482;

/// φ² - Phi squared
pub const PHI_SQ: f64 = 2.6180339887498948482;

/// 1/φ² - Inverse phi squared
pub const INV_PHI_SQ: f64 = 0.3819660112501051518;

/// 1/φ³ - Inverse phi cubed
pub const INV_PHI_CU: f64 = 0.2360679774997896964;

/// √2 - Square root of 2
pub const SQRT2: f64 = 1.4142135623730950488;

/// 1/√2 - Inverse square root of 2
pub const INV_SQRT2: f64 = 0.7071067811865475244;

/// ρ (rho) - Plastic constant
pub const PLASTIC: f64 = 1.3247179572447458000;

/// 1/ρ - Inverse plastic constant
pub const INV_PLASTIC: f64 = 0.7548776662466927600;

/// π - Pi
pub const PI: f64 = std::f64::consts::PI;

/// 2π - Two pi
pub const TWO_PI: f64 = 2.0 * PI;

/// Ring buffer size (power of 2)
pub const RING_SIZE: usize = 131072;

/// Ring buffer mask for wrapping
pub const RING_MASK: usize = RING_SIZE - 1;

/// Maximum number of simultaneous grains
pub const MAX_GRAINS: usize = 32;

/// Tiny value to avoid denormals
pub const TINY: f64 = 1.0e-30;

/// Default amplitude normalization
pub const AMP_NORM: f64 = 0.55;

/// Output drive for soft clipping
pub const OUT_DRIVE: f64 = 1.2;

/// Maximum event rate in Hz
pub const MAX_EVENT_RATE_HZ: f64 = 240.0;

/// Minimum base grain length in ms
pub const MIN_BASE_LENGTH_MS: f64 = 1.0;

/// Minimum grain duration in samples
pub const MIN_GRAIN_SAMPLES: f64 = 8.0;

/// Default sample rate
pub const DEFAULT_SR: f64 = 44100.0;

// φ-based coefficients for hemisphere coupling (corrected from C++ magic numbers)
/// Minimum hemisphere weight (should be 1/φ² ≈ 0.382)
pub const HEMI_WEIGHT_MIN: f64 = INV_PHI_SQ;

/// Maximum hemisphere weight (should be 1/φ ≈ 0.618)  
pub const HEMI_WEIGHT_MAX: f64 = INV_PHI;

// Hawkes process φ-corrected parameters
/// Hawkes boost factor (1/φ instead of 0.7)
pub const HAWKES_BOOST: f64 = INV_PHI;
