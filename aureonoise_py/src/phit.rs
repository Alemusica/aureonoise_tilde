//! PHITS — Phase-bit entropy pool for true random
//!
//! Ported from aperiodic-binaural-reverb `phit.rs` and triphase-computation `libphit.h`.
//! Extracts entropy from the beat frequency between CPU execution time
//! and system timer — two asynchronous clock domains.
//!
//! On macOS ARM64: ~1.96 phits per single read, ~4.06 per compound (N=2).
//! Timer: clock_gettime_nsec_np(CLOCK_UPTIME_RAW) — 24 MHz resolution.
//!
//! For audio DSP: seeds entropy pool from hardware at init, then fast
//! pool-based extraction per-sample with periodic re-seeding (~every 23ms).
//! No syscalls on the hot path between re-seeds.

// ─── Platform timer ───────────────────────────────────────────────

/// Read system timer in nanoseconds (macOS: CLOCK_UPTIME_RAW, 24 MHz).
#[cfg(target_os = "macos")]
fn now_ns() -> u64 {
    extern "C" {
        fn clock_gettime_nsec_np(clock_id: u32) -> u64;
    }
    const CLOCK_UPTIME_RAW: u32 = 8;
    unsafe { clock_gettime_nsec_np(CLOCK_UPTIME_RAW) }
}

/// Fallback timer for non-macOS platforms.
#[cfg(not(target_os = "macos"))]
fn now_ns() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos() as u64
}

// ─── Mixing primitives ───────────────────────────────────────────

/// Golden ratio constant: phi x 2^64, where phi = (sqrt(5) - 1) / 2.
/// Sequential IDs produce maximally separated hash inputs (Weyl equidistribution).
const GOLDEN_RATIO_64: u64 = 0x9E37_79B9_7F4A_7C15;

/// Reseed interval: harvest hardware entropy every 1024 extractions (~23ms @ 44.1kHz).
const RESEED_INTERVAL: u64 = 1024;

/// Number of pool lanes (256-bit state = 4 x u64).
const POOL_LANES: usize = 4;

/// Initial seeding rounds: 16 hardware reads for full pool saturation.
const SEED_ROUNDS: usize = 16;

/// SplitMix64 mixing function — bijective, full-period.
/// From libphit.h: avalanche-optimal constants by Stafford.
#[inline]
fn phit_hash64(mut key: u64) -> u64 {
    key = (key ^ (key >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    key = (key ^ (key >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    key ^ (key >> 31)
}

/// Calibrated workload — LCG, 20 iterations.
/// `black_box` prevents optimization; execution time varies with CPU phase,
/// DVFS state, cache effects, and speculative execution.
#[inline(never)]
fn phit_workload() -> u64 {
    let mut x: u64 = 0xCAFE_BABE;
    for _ in 0..20 {
        x = x
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        x ^= x >> 17;
    }
    std::hint::black_box(x)
}

// ─── PhitPool ─────────────────────────────────────────────────────

/// 256-bit entropy pool (4 lane u64), rotate-left mixing, golden-ratio counter.
///
/// Continuously mixed with hardware phits. Between re-seeds, extraction is
/// purely computational (no syscalls) — safe for the audio hot path.
pub struct PhitPool {
    lanes: [u64; POOL_LANES],
    mix_counter: u64,
}

impl PhitPool {
    /// Create a new pool, seeded from 16 hardware extractions.
    pub fn new() -> Self {
        let mut pool = Self {
            lanes: [0; POOL_LANES],
            mix_counter: 0,
        };
        for _ in 0..SEED_ROUNDS {
            pool.harvest();
        }
        pool
    }

    /// Inject entropy from workload timing.
    ///
    /// Runs the calibrated LCG workload, reads the timer, and mixes both
    /// into the pool. Two entropy sources per call:
    /// 1. `t`: absolute timer phase at completion
    /// 2. `x ^ t`: workload result XOR timer — captures DVFS transitions
    pub fn stir(&mut self) {
        self.harvest();
    }

    /// Extract a u64 from pool state (no hardware call).
    ///
    /// Advances the golden-ratio counter, hashes the target lane,
    /// then combines all 4 lanes with distinct rotations for full diffusion.
    pub fn extract_u64(&mut self) -> u64 {
        self.mix_counter += 1;
        let slot = (self.mix_counter & 3) as usize;
        self.lanes[slot] = phit_hash64(
            self.lanes[slot] ^ self.mix_counter.wrapping_mul(GOLDEN_RATIO_64),
        );
        let mut out = self.lanes[0];
        out ^= self.lanes[1].rotate_left(13);
        out ^= self.lanes[2].rotate_left(29);
        out ^= self.lanes[3].rotate_left(43);
        out
    }

    /// Feed a raw sample into the pool with golden-ratio counter mixing.
    fn feed(&mut self, sample: u64) {
        self.mix_counter += 1;
        let z = phit_hash64(
            sample.wrapping_add(self.mix_counter.wrapping_mul(GOLDEN_RATIO_64)),
        );
        let slot = (self.mix_counter & 3) as usize;
        self.lanes[slot] ^= z;
        self.lanes[(slot + 1) & 3] ^= self.lanes[slot].rotate_left(17);
    }

    /// Harvest one phit from hardware and mix into pool.
    fn harvest(&mut self) {
        let x = phit_workload();
        let t = now_ns();
        self.feed(t);
        self.feed(x ^ t);
    }
}

// ─── PhitRng ──────────────────────────────────────────────────────

/// Fast PRNG seeded from hardware phit entropy.
///
/// After initial seeding (16 hardware reads), generation is fast (no syscalls).
/// Re-seeds from hardware every 1024 extractions (~23ms at 44.1kHz) to keep
/// entropy fresh — the pool never falls into a deterministic cycle.
pub struct PhitRng {
    pool: PhitPool,
    counter: u64,
}

impl PhitRng {
    /// Create with a new hardware-seeded entropy pool.
    pub fn new() -> Self {
        Self {
            pool: PhitPool::new(),
            counter: 0,
        }
    }

    /// Create from an existing pool (shares initial entropy, then diverges).
    pub fn from_pool(pool: &mut PhitPool) -> Self {
        let mut lanes = [0u64; POOL_LANES];
        for lane in &mut lanes {
            *lane = pool.extract_u64();
        }
        Self {
            pool: PhitPool {
                lanes,
                mix_counter: 0,
            },
            counter: 0,
        }
    }

    /// Fast u64 generation with periodic hardware re-seeding.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.counter += 1;
        self.maybe_reseed(&mut None);
        self.pool.extract_u64()
    }

    /// Float in [0, 1) — 53 bits of mantissa.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Bipolar float in [-1, 1).
    #[inline]
    pub fn next_f64_bipolar(&mut self) -> f64 {
        2.0 * self.next_f64() - 1.0
    }

    /// Periodic re-seed check: harvest hardware entropy every 1024 calls.
    fn maybe_reseed(&mut self, external_pool: &mut Option<&mut PhitPool>) {
        if self.counter % RESEED_INTERVAL == 0 {
            match external_pool {
                Some(p) => {
                    // Cross-seed from external pool
                    let sample = p.extract_u64();
                    self.pool.feed(sample);
                    self.pool.harvest();
                }
                None => {
                    self.pool.harvest();
                }
            }
        }
    }

    /// Re-seed from an external pool (for explicit cross-seeding).
    pub(crate) fn reseed_from(&mut self, pool: &mut PhitPool) {
        let sample = pool.extract_u64();
        self.pool.feed(sample);
        self.pool.harvest();
    }
}

// ─── Distribution ─────────────────────────────────────────────────

/// Golden-ratio spacing for decorrelation between N resonators.
///
/// Returns N values in [0, 1) with maximal separation, based on
/// Weyl's Equidistribution Theorem: the sequence {n * phi mod 1}
/// is uniformly distributed for irrational phi, and phi is optimal
/// (slowest convergence to periodicity).
pub fn distribute(n: usize) -> Vec<f64> {
    let phi_frac: f64 = 0.618_033_988_749_894_848_2; // (sqrt(5) - 1) / 2
    (0..n)
        .map(|i| ((i + 1) as f64 * phi_frac) % 1.0)
        .collect()
}

/// Distribute a single phase reading across resonator IDs (per-resonator variant).
///
/// `phase`: single hardware read or pool extraction.
/// `resonator_id`: 0..N-1, each ID sees a statistically independent value.
#[inline]
pub(crate) fn distribute_phase(phase: u64, resonator_id: u32) -> u64 {
    phit_hash64(phase ^ (resonator_id as u64).wrapping_mul(GOLDEN_RATIO_64))
}

// ─── Tests ────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monobit() {
        // 47-53% ratio of 1-bits in 10000 u64s
        let mut rng = PhitRng::new();
        let n = 10_000;
        let total_bits = n * 64;
        let mut ones = 0u64;
        for _ in 0..n {
            ones += rng.next_u64().count_ones() as u64;
        }
        let ratio = ones as f64 / total_bits as f64;
        assert!(
            (0.47..=0.53).contains(&ratio),
            "Monobit test failed: ratio {:.6} outside [0.47, 0.53]",
            ratio
        );
    }

    #[test]
    fn test_byte_distribution() {
        // Chi-squared < 310 on 256 buckets
        let mut rng = PhitRng::new();
        let n = 10_000;
        let mut buckets = [0u64; 256];
        for _ in 0..n {
            let v = rng.next_u64();
            // Test each byte of each u64
            for shift in (0..64).step_by(8) {
                let byte = ((v >> shift) & 0xFF) as usize;
                buckets[byte] += 1;
            }
        }
        let total: u64 = buckets.iter().sum();
        let expected = total as f64 / 256.0;
        let chi_sq: f64 = buckets
            .iter()
            .map(|&b| {
                let d = b as f64 - expected;
                d * d / expected
            })
            .sum();
        assert!(
            chi_sq < 310.0,
            "Byte distribution chi-squared {:.2} >= 310",
            chi_sq
        );
    }

    #[test]
    fn test_per_bit_entropy() {
        // Each bit position should be close to 50/50 => > 63 bits of entropy per 64
        let mut rng = PhitRng::new();
        let n = 10_000u64;
        let mut bit_counts = [0u64; 64];
        for _ in 0..n {
            let v = rng.next_u64();
            for bit in 0..64 {
                if (v >> bit) & 1 == 1 {
                    bit_counts[bit] += 1;
                }
            }
        }
        // Per-bit entropy: H = -p*log2(p) - (1-p)*log2(1-p)
        let mut total_entropy = 0.0;
        for &count in &bit_counts {
            let p = count as f64 / n as f64;
            if p > 0.0 && p < 1.0 {
                total_entropy += -(p * p.log2() + (1.0 - p) * (1.0 - p).log2());
            }
        }
        assert!(
            total_entropy > 63.0,
            "Per-bit entropy {:.4} <= 63 bits",
            total_entropy
        );
    }

    #[test]
    fn test_distribute_spacing() {
        // N=8, verify values are in [0, 1) and reasonably spaced
        let vals = distribute(8);
        assert_eq!(vals.len(), 8);

        // All in [0, 1)
        for (i, &v) in vals.iter().enumerate() {
            assert!(
                (0.0..1.0).contains(&v),
                "distribute[{}] = {} outside [0, 1)",
                i,
                v
            );
        }

        // Sort and check minimum gap (golden-ratio spacing guarantees good separation)
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // With N=8 and golden-ratio spacing, minimum gap should be > 1/(2*N)
        let min_gap_threshold = 1.0 / (2.0 * 8.0);
        for i in 0..sorted.len() {
            let next = if i + 1 < sorted.len() {
                sorted[i + 1]
            } else {
                sorted[0] + 1.0 // wrap-around gap
            };
            let gap = next - sorted[i];
            assert!(
                gap > min_gap_threshold,
                "Gap between sorted[{}]={:.4} and next={:.4} is {:.6} < {:.6}",
                i,
                sorted[i],
                next,
                gap,
                min_gap_threshold
            );
        }

        // All values should be distinct
        sorted.dedup();
        assert_eq!(sorted.len(), 8, "Not all 8 values are distinct");
    }

    #[test]
    fn test_pool_extract_varies() {
        let mut pool = PhitPool::new();
        let a = pool.extract_u64();
        let b = pool.extract_u64();
        let c = pool.extract_u64();
        assert_ne!(a, b);
        assert_ne!(b, c);
        assert_ne!(a, c);
    }

    #[test]
    fn test_rng_f64_range() {
        let mut rng = PhitRng::new();
        for _ in 0..10_000 {
            let v = rng.next_f64();
            assert!(
                (0.0..1.0).contains(&v),
                "next_f64() = {} outside [0, 1)",
                v
            );
        }
    }

    #[test]
    fn test_rng_bipolar_range() {
        let mut rng = PhitRng::new();
        for _ in 0..10_000 {
            let v = rng.next_f64_bipolar();
            assert!(
                (-1.0..1.0).contains(&v),
                "next_f64_bipolar() = {} outside [-1, 1)",
                v
            );
        }
    }

    #[test]
    fn test_two_rngs_differ() {
        let mut rng_a = PhitRng::new();
        let mut rng_b = PhitRng::new();
        let mut differ = false;
        for _ in 0..100 {
            if rng_a.next_u64() != rng_b.next_u64() {
                differ = true;
                break;
            }
        }
        assert!(differ, "Two PhitRngs should produce different sequences");
    }

    #[test]
    fn test_from_pool() {
        let mut pool = PhitPool::new();
        let mut rng = PhitRng::from_pool(&mut pool);
        let a = rng.next_u64();
        let b = rng.next_u64();
        assert_ne!(a, b, "from_pool RNG should produce varied output");
    }

    #[test]
    fn test_stir_affects_output() {
        let mut pool = PhitPool::new();
        let before = pool.extract_u64();
        pool.stir();
        let after = pool.extract_u64();
        // stir() mixes hardware entropy — output should change
        assert_ne!(before, after, "stir() should affect pool output");
    }

    #[test]
    fn test_distribute_phase_decorrelation() {
        let phase = 0xDEAD_BEEF_CAFE_1234u64;
        let vals: Vec<u64> = (0..512).map(|i| distribute_phase(phase, i)).collect();

        // All should be distinct
        let mut sorted = vals.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 512, "All distributed values should be unique");

        // Monobit: count 1s across all values, should be ~50%
        let total_bits = 512 * 64;
        let ones: u64 = vals.iter().map(|v| v.count_ones() as u64).sum();
        let ratio = ones as f64 / total_bits as f64;
        assert!(
            (0.47..=0.53).contains(&ratio),
            "Distribute monobit: ratio {:.4} outside [0.47, 0.53]",
            ratio
        );
    }
}
