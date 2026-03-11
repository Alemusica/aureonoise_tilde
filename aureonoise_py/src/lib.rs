//! aureonoise - Main DSP engine
//! Granular noise/glitch texture generator

use pyo3::prelude::*;
use numpy::PyArray1;

mod constants;
mod rng;
mod weyl;
mod math;
mod noise;
mod ring;
mod stoch;
mod envelope;
mod grain;
mod burst;
mod phi_model;
mod phit;
mod sr;
mod external;
mod modal;
mod dialogue;
mod binaural;
mod isochronic;
mod tinnitus;
mod dvf;
mod room;
mod polyrhythm;
mod phi_lattice;

pub use constants::*;
pub use rng::Rng;
pub use weyl::Weyl;
pub use math::*;
pub use noise::{NoiseColor, NoiseColorState, NoiseGen, NoiseMode, GrainKind};
pub use ring::RingBuffer;
pub use stoch::{OrnsteinUhlenbeck, Lattice, Hawkes};
pub use envelope::{Envelope, EnvelopeShape};
pub use grain::{Grain, GrainPool};
pub use burst::{BurstEngine, BurstResult};
pub use phi_model::{PhiModel, GeometryConfig, Geometry, build_geometry, compute_head_result,
    design_pinna_response, PinnaTuning, compute_distance_response, compute_air_absorption};
pub use external::{ExternalProcessor, ExternalConfig};
pub use dialogue::{DialogueSystem, DialogueParams, PhiPan, BilateralOscillator};
pub use modal::{ModalEngine, ModalPreset};
pub use binaural::BinauralBeat;
pub use isochronic::IsochronicTone;
pub use tinnitus::TinnitusNotch;
pub use noise::SpectralTilt;
pub use dvf::DvfFilter;
pub use room::RoomReverb;
pub use polyrhythm::PolyrhythmClock;
pub use phi_lattice::PhiLattice;
pub use sr::StochasticResonance;

/// aureonoise DSP engine parameters
#[pyclass]
#[derive(Clone)]
pub struct Params {
    // Timing
    #[pyo3(get, set)]
    pub rate: f64,
    #[pyo3(get, set)]
    pub baselen_ms: f64,
    #[pyo3(get, set)]
    pub len_phi: f64,
    
    // Spatial
    #[pyo3(get, set)]
    pub width: f64,
    #[pyo3(get, set)]
    pub itd_us: f64,
    #[pyo3(get, set)]
    pub ild_db: f64,
    #[pyo3(get, set)]
    pub hemis_coupling: f64,
    #[pyo3(get, set)]
    pub spat_min_deg: f64,
    #[pyo3(get, set)]
    pub spat_min_ms: f64,
    #[pyo3(get, set)]
    pub spat_ipd: f64,
    #[pyo3(get, set)]
    pub spat_shadow: f64,
    
    // Envelope
    #[pyo3(get, set)]
    pub env_attack: f64,
    #[pyo3(get, set)]
    pub env_decay: f64,
    #[pyo3(get, set)]
    pub env_sustain: f64,
    #[pyo3(get, set)]
    pub env_release: f64,
    /// Envelope curve shape: 0 = linear (default), 1 = Hann (raised cosine attack/release)
    #[pyo3(get, set)]
    pub envelope_shape: i32,

    // Timbre
    #[pyo3(get, set)]
    pub noise_color: i32,
    #[pyo3(get, set)]
    pub color_amt: f64,
    #[pyo3(get, set)]
    pub glitch_mix: f64,
    #[pyo3(get, set)]
    pub srcrush_amt: f64,
    #[pyo3(get, set)]
    pub bitcrush_amt: f64,
    
    // Stochastic
    #[pyo3(get, set)]
    pub thermo: bool,
    #[pyo3(get, set)]
    pub lattice: bool,
    #[pyo3(get, set)]
    pub burst: bool,
    #[pyo3(get, set)]
    pub temperature: f64,
    #[pyo3(get, set)]
    pub lat_rate: f64,
    #[pyo3(get, set)]
    pub lat_eps: f64,
    #[pyo3(get, set)]
    pub lat_gamma: f64,
    #[pyo3(get, set)]
    pub lat_sigma: f64,
    #[pyo3(get, set)]
    pub burst_floor: f64,
    #[pyo3(get, set)]
    pub burst_phi_mix: f64,

    // Stochastic resonance (Collins 1995)
    #[pyo3(get, set)]
    pub sr_on: bool,

    // Spatial — externalisation
    #[pyo3(get, set)]
    pub externalization: f64,

    // Dialogue (interhemispheric coherence)
    #[pyo3(get, set)]
    pub dialogue_on: bool,
    #[pyo3(get, set)]
    pub dialogue_strength: f64,
    #[pyo3(get, set)]
    pub dialogue_memory: f64,
    #[pyo3(get, set)]
    pub dialogue_phi_mix: f64,

    // Phi-Pan + Bilateral
    #[pyo3(get, set)]
    pub phi_pan: bool,
    #[pyo3(get, set)]
    pub bilateral_on: bool,
    #[pyo3(get, set)]
    pub bilateral_rate: f64,
    #[pyo3(get, set)]
    pub bilateral_amount: f64,
    /// Theta-gamma nesting: lock grain rate as integer multiple of bilateral rate.
    /// Lisman-Jensen 2013: 4-8 gamma cycles per theta cycle.
    #[pyo3(get, set)]
    pub bilateral_nesting: bool,

    // Noise (extended modes 0-5)
    #[pyo3(get, set)]
    pub noise_mode: i32,
    #[pyo3(get, set)]
    pub aureo_decay: f64,
    #[pyo3(get, set)]
    pub aureo_stride: f64,
    #[pyo3(get, set)]
    pub aureo_harmonics: i32,
    #[pyo3(get, set)]
    pub quantum_detail: f64,
    #[pyo3(get, set)]
    pub quantum_base: f64,
    #[pyo3(get, set)]
    pub velvet_density: f64,

    // Modal resonator
    #[pyo3(get, set)]
    pub modal_on: bool,
    #[pyo3(get, set)]
    pub modal_mix: f64,
    #[pyo3(get, set)]
    pub modal_decay: f64,
    #[pyo3(get, set)]
    pub modal_preset: i32,
    #[pyo3(get, set)]
    pub modal_mirror: f64,
    #[pyo3(get, set)]
    pub modal_feedback: f64,
    #[pyo3(get, set)]
    pub modal_contralateral: f64,  // 0.0-1.0, strength of spatial mirror

    // Coherence feedback (BAC-inspired closed-loop)
    /// Enable coherence feedback loop (opt-in for therapeutic presets)
    #[pyo3(get, set)]
    pub feedback_on: bool,
    /// Temperature ramp duration in seconds (0 = instant, >0 = linear ramp)
    #[pyo3(get, set)]
    pub temp_ramp_sec: f64,

    // Binaural beat
    #[pyo3(get, set)]
    pub binaural_on: bool,
    #[pyo3(get, set)]
    pub binaural_carrier_hz: f64,
    #[pyo3(get, set)]
    pub binaural_beat_hz: f64,
    #[pyo3(get, set)]
    pub binaural_level: f64,

    // Isochronic tone
    #[pyo3(get, set)]
    pub isochronic_on: bool,
    #[pyo3(get, set)]
    pub isochronic_carrier_hz: f64,
    #[pyo3(get, set)]
    pub isochronic_rate_hz: f64,
    #[pyo3(get, set)]
    pub isochronic_duty: f64,
    #[pyo3(get, set)]
    pub isochronic_level: f64,

    // Continuous spectral slope (replaces discrete noise_color for tilt)
    #[pyo3(get, set)]
    pub noise_slope: f64,

    // Tinnitus notch (0=off, >0=center freq Hz)
    #[pyo3(get, set)]
    pub tinnitus_notch_hz: f64,
    #[pyo3(get, set)]
    pub tinnitus_notch_q: f64,

    // Phi model (expose for GUI)
    #[pyo3(get, set)]
    pub phi_distance: f64,
    #[pyo3(get, set)]
    pub phi_elev: f64,
    /// Pinna amplitude modulation amount (0-1). 0 = bypass, 1 = full pinna response.
    #[pyo3(get, set)]
    pub spat_pinna: f64,
    /// Distance attenuation + LP rolloff amount (0-1). 0 = bypass, 1 = full distance model.
    #[pyo3(get, set)]
    pub spat_distance: f64,

    // Polyrhythm clock (T5.5)
    #[pyo3(get, set)]
    pub polyrhythm_on: bool,
    #[pyo3(get, set)]
    pub polyrhythm_p: u32,
    #[pyo3(get, set)]
    pub polyrhythm_q: u32,
    #[pyo3(get, set)]
    pub polyrhythm_rate: f64,
    #[pyo3(get, set)]
    pub polyrhythm_amount: f64,

    // Room reverb (T6.2)
    #[pyo3(get, set)]
    pub room_mix: f64,

    // Coherence-driven spatial morphing (T6.3)
    #[pyo3(get, set)]
    pub coherence_spatial: bool,

    // Phi lattice (universal phi ratio oracle)
    #[pyo3(get, set)]
    pub phi_lattice_on: bool,
    /// Per-grain personality amount (0=all same, 1=max variation)
    #[pyo3(get, set)]
    pub phi_personality: f64,
    /// Timing snap to phi grid (0=free, 1=strict)
    #[pyo3(get, set)]
    pub phi_timing_strength: f64,
    /// Spatial snap to phi lattice (0=free, 1=strict)
    #[pyo3(get, set)]
    pub phi_spatial_strength: f64,

    // System
    #[pyo3(get, set)]
    pub seed: u64,
}

#[pymethods]
impl Params {
    #[new]
    pub fn new() -> Self {
        Self::default()
    }
}

impl Default for Params {
    fn default() -> Self {
        Self {
            // Timing
            rate: 8.0,
            baselen_ms: 120.0,
            len_phi: 0.8,
            
            // Spatial
            width: 1.0,
            itd_us: 600.0,
            ild_db: 6.0,
            hemis_coupling: 0.6,
            spat_min_deg: 12.0,
            spat_min_ms: 35.0,
            spat_ipd: 0.6,
            spat_shadow: 0.7,
            
            // Envelope
            env_attack: 0.18,
            env_decay: 0.28,
            env_sustain: 0.55,
            env_release: 0.30,
            envelope_shape: 0,  // linear (default)

            // Timbre
            noise_color: 1, // Pink
            color_amt: 0.65,
            glitch_mix: 0.5,
            srcrush_amt: 0.2,
            bitcrush_amt: 0.15,
            
            // Stochastic
            thermo: true,
            lattice: true,
            burst: true,
            temperature: 0.22,
            lat_rate: 250.0,
            lat_eps: INV_PHI_CU,
            lat_gamma: PHI,
            lat_sigma: 0.06,
            burst_floor: 0.35,
            burst_phi_mix: 0.6,

            // Stochastic resonance
            sr_on: false,

            // Spatial — externalisation
            externalization: 0.0,

            // Dialogue
            dialogue_on: true,
            dialogue_strength: 0.6,
            dialogue_memory: 0.5,
            dialogue_phi_mix: 0.75,

            // Phi-Pan + Bilateral
            phi_pan: false,
            bilateral_on: false,
            bilateral_rate: 1.0,
            bilateral_amount: 0.8,
            bilateral_nesting: false,

            // Noise (extended modes)
            noise_mode: 1,  // Pink
            aureo_decay: 0.3,
            aureo_stride: 1.0,
            aureo_harmonics: 12,
            quantum_detail: 0.7,
            quantum_base: 220.0,
            velvet_density: 2000.0,

            // Modal
            modal_on: false,
            modal_mix: 0.3,
            modal_decay: 0.5,
            modal_preset: 1, // Wood
            modal_mirror: 0.3,
            modal_feedback: 0.1,
            modal_contralateral: 0.0,

            // Feedback
            feedback_on: false,
            temp_ramp_sec: 0.0,

            // Binaural beat
            binaural_on: false,
            binaural_carrier_hz: 250.0,
            binaural_beat_hz: 6.0,
            binaural_level: 0.08,

            // Isochronic
            isochronic_on: false,
            isochronic_carrier_hz: 165.0,
            isochronic_rate_hz: 10.0,
            isochronic_duty: 0.5,
            isochronic_level: 0.10,

            // Spectral slope
            noise_slope: -1.0, // pink default (backward compat)

            // Tinnitus notch
            tinnitus_notch_hz: 0.0, // off
            tinnitus_notch_q: 6.0,

            // Phi model
            phi_distance: 1.5,
            phi_elev: 0.0,
            spat_pinna: 0.0,
            spat_distance: 0.0,

            // Polyrhythm clock
            polyrhythm_on: false,
            polyrhythm_p: 3,
            polyrhythm_q: 2,
            polyrhythm_rate: 0.5,
            polyrhythm_amount: 0.5,

            // Room reverb
            room_mix: 0.0,

            // Coherence spatial morphing
            coherence_spatial: false,

            // Phi lattice
            phi_lattice_on: false,
            phi_personality: 0.5,
            phi_timing_strength: 0.5,
            phi_spatial_strength: 0.5,

            // System
            seed: 20251010,
        }
    }
}

/// aureonoise DSP engine
#[pyclass]
pub struct Engine {
    params: Params,
    sr: f64,

    // State
    rng: Rng,
    w_phi: Weyl,
    w_s2: Weyl,
    w_pl: Weyl,
    noise_gen: NoiseGen,
    ring: RingBuffer,
    grains: GrainPool,
    envelope: Envelope,

    // Burst position modulation
    burst_engine: BurstEngine,
    burst_centroid: f64,      // EMA of burst-weighted grain pan positions
    burst_centroid_alpha: f64, // EMA smoothing factor (~25ms time constant)

    // External externalisation (block-level cross-channel feedback delay)
    external_proc: ExternalProcessor,

    // Dialogue (interhemispheric coherence)
    dialogue: DialogueSystem,
    phi_pan_proc: PhiPan,
    bilateral: BilateralOscillator,

    // Modal resonator
    modal_engine: ModalEngine,

    // Phi head geometry (for ITD computation)
    phi_geom: Geometry,
    phi_itd_max: f64, // max ITD samples at pan=1.0 for normalization

    // New DSP modules (Sprint 4)
    binaural: BinauralBeat,
    isochronic: IsochronicTone,
    tinnitus: TinnitusNotch,
    spectral_tilt: SpectralTilt,

    // Sprint 6 modules
    dvf: DvfFilter,
    room: RoomReverb,
    polyrhythm: PolyrhythmClock,

    // Stochastic resonance (Collins 1995)
    stoch_res: StochasticResonance,

    // Phi lattice (universal phi ratio oracle)
    phi_lattice: phi_lattice::PhiLattice,

    // Hardware-seeded entropy (phit)
    phit_rng: phit::PhitRng,

    // Stochastic
    ou_pan: OrnsteinUhlenbeck,
    ou_itd: OrnsteinUhlenbeck,
    ou_amp: OrnsteinUhlenbeck,
    ou_rate: OrnsteinUhlenbeck,
    lattice: Lattice,
    hawkes: Hawkes,
    lat_phase: f64,
    lat_last_v: f64,

    // Scheduling
    samples_to_next: i32,
    gap_elapsed: i32,
    sample_counter: u64,

    // 0.1 Hz macro-modulation OU (heart-brain coherence rhythm)
    ou_macro: OrnsteinUhlenbeck,

    // Coherence feedback loop (BAC-inspired)
    coh_slow: f64,           // slow EMA of coherence (~5s tau)
    effective_temp: f64,     // feedback-modulated temperature
    effective_bilateral_rate: f64, // feedback-modulated bilateral rate

    // Temperature ramp
    temp_ramp_samples: f64,  // total ramp duration in samples (0 = instant)
    temp_ramp_elapsed: f64,  // samples elapsed in current ramp
    temp_start: f64,         // starting temperature
    temp_target: f64,        // target temperature

    // Weighted-average pan per block (for DVF near-field)
    block_pan_sum: f64,
    block_pan_weight: f64,

    // Contralateral modal spatial state
    contra_delay: [f64; 64],   // mono delay ring for ITD
    contra_delay_pos: usize,
    contra_shadow_z_l: f64,    // head shadow LP state (left ear)
    contra_shadow_z_r: f64,    // head shadow LP state (right ear)

    // Body resonance safety notch filters (P0 safety)
    // Each filter: [b0, b1, b2, a1, a2] coefficients, [z1, z2] state per channel
    safety_notch_chest_coeff: [f64; 5],   // 6.5 Hz center (5-8 Hz chest cavity)
    safety_notch_chest_zl: [f64; 2],      // left channel state
    safety_notch_chest_zr: [f64; 2],      // right channel state
    safety_notch_eye_coeff: [f64; 5],     // 19 Hz center (eyeball resonance)
    safety_notch_eye_zl: [f64; 2],        // left channel state
    safety_notch_eye_zr: [f64; 2],        // right channel state

    // Previous grain state (for hemisphere coupling)
    prev_pan: f64,
    prev_itd: f64,
    prev_ild: f64,
    last_gap_samples: f64,
    last_dur_samples: f64,
}

#[pymethods]
impl Engine {
    #[new]
    #[pyo3(signature = (sample_rate = 44100.0))]
    pub fn new(sample_rate: f64) -> Self {
        let sr = if sample_rate > 0.0 { sample_rate } else { 44100.0 };
        let params = Params::default();
        
        let mut rng = Rng::new(params.seed);
        let w_phi = Weyl::phi(rng.uni01());
        let w_s2 = Weyl::phi_sq(rng.uni01());
        let w_pl = Weyl::phi_cu(rng.uni01());
        
        let mut noise_gen = NoiseGen::new(NoiseMode::Pink);
        noise_gen.classic.set_amount(params.color_amt);

        let mut modal_engine = ModalEngine::new(sr);
        modal_engine.set_preset(ModalPreset::from_i32(params.modal_preset));
        modal_engine.set_active(params.modal_on);

        let phi_geom = build_geometry(&GeometryConfig::default());
        let phi_itd_max = compute_head_result(&phi_geom, sr, 1.0, 0.0, 1.0)
            .itd_samples.abs().max(1e-9);

        let init_temp = params.temperature;
        let init_bilateral_rate = params.bilateral_rate;

        // Body resonance safety notch filters (P0)
        // Chest cavity: 6.5 Hz center, Q=2.5 → bandwidth ~2.6 Hz (covers 5.2-7.8 Hz, avoids 3 Hz delta)
        // Eyeball: 19 Hz center, Q=4.0 → bandwidth ~4.75 Hz (covers 16.6-21.4 Hz)
        let safety_notch_chest_coeff = compute_notch_coeffs(6.5, sr, 2.5);
        let safety_notch_eye_coeff = compute_notch_coeffs(19.0, sr, 4.0);

        Self {
            params,
            sr,
            rng,
            w_phi,
            w_s2,
            w_pl,
            noise_gen,
            ring: RingBuffer::new(),
            grains: GrainPool::new(),
            envelope: Envelope::new(),
            burst_engine: BurstEngine {
                enabled: true,
                floor: 0.35,
                phi_mix: 0.6,
            },
            burst_centroid: 0.0,
            burst_centroid_alpha: 1.0 / (0.025 * sr),  // 25ms tau at sample rate
            external_proc: ExternalProcessor::new(),
            dialogue: DialogueSystem::new(),
            phi_pan_proc: PhiPan::new(),
            bilateral: BilateralOscillator::new(),
            modal_engine,
            phi_geom,
            phi_itd_max,
            binaural: BinauralBeat::new(),
            isochronic: IsochronicTone::new(),
            tinnitus: TinnitusNotch::new(sr),
            spectral_tilt: SpectralTilt::new(),
            dvf: DvfFilter::new(),
            room: RoomReverb::new(sr),
            polyrhythm: PolyrhythmClock::new(),
            stoch_res: StochasticResonance::new(),
            phi_lattice: phi_lattice::PhiLattice::new(),
            phit_rng: phit::PhitRng::new(),
            ou_pan: OrnsteinUhlenbeck::new(0.60, 0.0),
            ou_itd: OrnsteinUhlenbeck::new(0.40, 0.0),
            ou_amp: OrnsteinUhlenbeck::new(0.80, 0.0),
            ou_rate: OrnsteinUhlenbeck::new(1.20, 0.0),
            // 0.1 Hz macro OU: tau = 1/(2*pi*0.1) ≈ 1.59s, sigma = 0.25
            ou_macro: OrnsteinUhlenbeck::new(1.59, 0.25),
            lattice: Lattice::new(8, 8, 4),
            hawkes: Hawkes::new(None, None),
            lat_phase: 0.0,
            lat_last_v: 0.0,
            samples_to_next: (sr * 0.05) as i32,
            gap_elapsed: 0,
            sample_counter: 0,
            coh_slow: 0.5,
            effective_temp: init_temp,
            effective_bilateral_rate: init_bilateral_rate,
            temp_ramp_samples: 0.0,
            temp_ramp_elapsed: 0.0,
            temp_start: init_temp,
            temp_target: init_temp,
            block_pan_sum: 0.0,
            block_pan_weight: 0.0,
            contra_delay: [0.0; 64],
            contra_delay_pos: 0,
            contra_shadow_z_l: 0.0,
            contra_shadow_z_r: 0.0,
            safety_notch_chest_coeff,
            safety_notch_chest_zl: [0.0; 2],
            safety_notch_chest_zr: [0.0; 2],
            safety_notch_eye_coeff,
            safety_notch_eye_zl: [0.0; 2],
            safety_notch_eye_zr: [0.0; 2],
            prev_pan: 0.0,
            prev_itd: 0.0,
            prev_ild: 0.0,
            last_gap_samples: 0.0,
            last_dur_samples: 0.0,
        }
    }

    /// Set parameters
    pub fn set_params(&mut self, params: Params) {
        self.params = params;

        // Noise — use noise_mode (0-5) with fallback to noise_color (0-2)
        let mode = if self.params.noise_mode >= 0 && self.params.noise_mode <= 5 {
            self.params.noise_mode
        } else {
            self.params.noise_color.clamp(0, 2)
        };
        let nm = match mode {
            0 => NoiseMode::White,
            2 => NoiseMode::Brown,
            3 => NoiseMode::Aureo,
            4 => NoiseMode::Quantum,
            5 => NoiseMode::Velvet,
            _ => NoiseMode::Pink,
        };
        self.noise_gen.set_mode(nm);
        self.noise_gen.classic.set_amount(self.params.color_amt);
        self.noise_gen.set_aureo_decay(self.params.aureo_decay);
        self.noise_gen.set_aureo_stride(self.params.aureo_stride);
        self.noise_gen.set_aureo_harmonics(self.params.aureo_harmonics);
        self.noise_gen.set_quantum_detail(self.params.quantum_detail);
        self.noise_gen.set_quantum_base(self.params.quantum_base);
        self.noise_gen.set_velvet_density(self.params.velvet_density);

        // Lattice
        self.lattice.eps = self.params.lat_eps;
        self.lattice.gamma = self.params.lat_gamma;
        self.lattice.sigma = self.params.lat_sigma;

        // Burst
        self.burst_engine.enabled = self.params.burst;
        self.burst_engine.floor = self.params.burst_floor;
        self.burst_engine.phi_mix = self.params.burst_phi_mix;

        // Scale Hawkes base rate with grain rate so burst clustering tracks frequency.
        // Default rate=8.0 maps to default base=4*PHI. Clamped to [0.5, 50] × PHI.
        let rate_scale = clamp(self.params.rate / 8.0, 0.1, 12.0);
        self.hawkes.base = 4.0 * PHI * rate_scale;
        // Decay rate scales inversely: faster events = faster decay to avoid pile-up
        self.hawkes.beta = (30.0 / PHI) * rate_scale.sqrt();

        // Bilateral oscillator
        self.bilateral.set_rate(self.params.bilateral_rate);
        self.bilateral.set_amount(self.params.bilateral_amount);

        // Temperature ramp
        if self.params.temp_ramp_sec > 0.0 {
            self.temp_ramp_samples = self.params.temp_ramp_sec * self.sr;
            self.temp_ramp_elapsed = 0.0;
            self.temp_start = self.effective_temp;
            self.temp_target = self.params.temperature;
        } else {
            self.effective_temp = self.params.temperature;
        }
        self.effective_bilateral_rate = self.params.bilateral_rate;

        // DVF near-field (T6.1) — distance from phi_distance param
        self.dvf.update_params(self.params.phi_distance, self.sr);

        // Room reverb (T6.2)
        self.room.set_mix(self.params.room_mix);

        // Polyrhythm clock (T5.5)
        self.polyrhythm.p = self.params.polyrhythm_p;
        self.polyrhythm.q = self.params.polyrhythm_q;
        self.polyrhythm.base_rate = self.params.polyrhythm_rate;
        self.polyrhythm.active = self.params.polyrhythm_on;

        // Modal
        self.modal_engine.set_preset(ModalPreset::from_i32(self.params.modal_preset));
        self.modal_engine.set_active(self.params.modal_on);
        self.modal_engine.set_mix(self.params.modal_mix);
        self.modal_engine.set_decay_scale(self.params.modal_decay);
        self.modal_engine.set_mirror(self.params.modal_mirror);
        self.modal_engine.set_feedback(self.params.modal_feedback);
        self.modal_engine.contralateral = self.params.modal_contralateral;

        // Stochastic resonance
        self.stoch_res.set_active(self.params.sr_on);
    }

    /// Get current parameters
    pub fn get_params(&self) -> Params {
        self.params.clone()
    }

    /// Reset the engine
    pub fn reset(&mut self) {
        self.rng.seed(self.params.seed);
        self.w_phi = Weyl::phi(self.rng.uni01());
        self.w_s2 = Weyl::phi_sq(self.rng.uni01());
        self.w_pl = Weyl::phi_cu(self.rng.uni01());
        self.noise_gen.reset();
        self.ring.clear();
        self.grains.reset_all();
        self.lattice.reset();
        self.hawkes.reset();
        self.burst_engine.enabled = self.params.burst;
        self.burst_engine.floor = self.params.burst_floor;
        self.burst_engine.phi_mix = self.params.burst_phi_mix;
        self.burst_centroid = 0.0;
        self.external_proc.reset();
        self.dialogue.reset();
        self.phi_pan_proc.reset();
        self.bilateral.reset();
        self.modal_engine.reset();
        self.binaural.reset();
        self.isochronic.reset();
        self.tinnitus = TinnitusNotch::new(self.sr);
        self.spectral_tilt = SpectralTilt::new();
        self.dvf.reset();
        self.room.reset();
        self.polyrhythm.reset();
        self.stoch_res.reset();
        self.phi_lattice = phi_lattice::PhiLattice::new();
        self.ou_macro.reset();
        self.lat_phase = 0.0;
        self.lat_last_v = 0.0;
        self.samples_to_next = (self.sr * 0.05) as i32;
        self.gap_elapsed = 0;
        self.sample_counter = 0;
        self.coh_slow = 0.5;
        self.effective_temp = self.params.temperature;
        self.effective_bilateral_rate = self.params.bilateral_rate;
        self.temp_ramp_samples = self.params.temp_ramp_sec * self.sr;
        self.temp_ramp_elapsed = 0.0;
        self.temp_start = self.params.temperature;
        self.temp_target = self.params.temperature;
        self.block_pan_sum = 0.0;
        self.block_pan_weight = 0.0;
        self.contra_delay = [0.0; 64];
        self.contra_delay_pos = 0;
        self.contra_shadow_z_l = 0.0;
        self.contra_shadow_z_r = 0.0;
        // Recompute safety notch coefficients (SR may have changed)
        self.safety_notch_chest_coeff = compute_notch_coeffs(6.5, self.sr, 2.5);
        self.safety_notch_chest_zl = [0.0; 2];
        self.safety_notch_chest_zr = [0.0; 2];
        self.safety_notch_eye_coeff = compute_notch_coeffs(19.0, self.sr, 4.0);
        self.safety_notch_eye_zl = [0.0; 2];
        self.safety_notch_eye_zr = [0.0; 2];
        self.prev_pan = 0.0;
        self.prev_itd = 0.0;
        self.prev_ild = 0.0;
        self.last_gap_samples = 0.0;
        self.last_dur_samples = 0.0;
    }

    /// Get current coherence level from the dialogue system (for GUI)
    pub fn coherence(&self) -> f64 {
        self.dialogue.coherence()
    }

    /// Get handshake count from dialogue system
    pub fn handshake_count(&self) -> u64 {
        self.dialogue.handshake_count()
    }

    /// Get handshake ratio (handshakes / utterances)
    pub fn handshake_ratio(&self) -> f64 {
        self.dialogue.handshake_ratio()
    }

    /// Get mean coherence across all utterances
    pub fn coherence_mean(&self) -> f64 {
        self.dialogue.coherence_mean()
    }

    /// Get Phase Locking Value — rhythmicity of handshakes [0, 1].
    /// 1.0 = perfectly periodic handshakes, 0.0 = random timing.
    pub fn handshake_plv(&self) -> f64 {
        self.dialogue.handshake_plv()
    }

    /// Get current stochastic resonance noise gain (for GUI display).
    /// Returns 1.0 when SR is inactive.
    pub fn sr_noise_gain(&self) -> f64 {
        self.stoch_res.noise_gain()
    }

    /// Process a block of samples, returns (left, right) arrays
    pub fn process<'py>(
        &mut self,
        py: Python<'py>,
        num_samples: usize,
    ) -> (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>) {
        let mut out_l = vec![0.0f64; num_samples];
        let mut out_r = vec![0.0f64; num_samples];
        
        self.process_block(&mut out_l, &mut out_r);
        
        (
            PyArray1::from_slice_bound(py, &out_l),
            PyArray1::from_slice_bound(py, &out_r),
        )
    }

    /// Get sample rate
    pub fn sample_rate(&self) -> f64 {
        self.sr
    }

    /// Set sample rate
    pub fn set_sample_rate(&mut self, sr: f64) {
        self.sr = if sr > 0.0 { sr } else { 44100.0 };
        self.modal_engine.set_sr(self.sr);
        // Recompute safety notch coefficients for new SR
        self.safety_notch_chest_coeff = compute_notch_coeffs(6.5, self.sr, 2.5);
        self.safety_notch_eye_coeff = compute_notch_coeffs(19.0, self.sr, 4.0);
    }

    /// Returns true when the current isochronic rate is in the 8-25 Hz
    /// seizure risk range (auditory driving, analogous to photic driving).
    /// The Python layer should show a safety warning when this returns true.
    pub fn isochronic_seizure_risk(&self) -> bool {
        self.params.isochronic_on
            && IsochronicTone::is_seizure_risk_range(self.params.isochronic_rate_hz)
    }

    /// Returns true when the current isochronic rate is in a body resonance
    /// range: 5-8 Hz (thoracic cavity) or 18-20 Hz (ocular globe).
    /// The Python layer should show a safety warning when this returns true.
    pub fn isochronic_body_resonance_risk(&self) -> bool {
        self.params.isochronic_on
            && IsochronicTone::is_body_resonance_range(self.params.isochronic_rate_hz)
    }
}

impl Engine {
    /// Internal block processing
    pub fn process_block(&mut self, out_l: &mut [f64], out_r: &mut [f64]) {
        let num_samples = out_l.len().min(out_r.len());

        // --- Temperature ramp (T4.2) ---
        // Ramp computes baseline; feedback modulates on top. No overwrite conflict.
        let ramped_temp = if self.temp_ramp_samples > 0.0 && self.temp_ramp_elapsed < self.temp_ramp_samples {
            self.temp_ramp_elapsed += num_samples as f64;
            let t = clamp01(self.temp_ramp_elapsed / self.temp_ramp_samples);
            self.temp_start + t * (self.temp_target - self.temp_start)
        } else {
            self.effective_temp
        };

        // --- Coherence feedback loop (T4.1, BAC-inspired) ---
        if self.params.feedback_on {
            let coh = self.dialogue.coherence();
            // Normalize: coherence typically in [0.6, 1.8], map to [0, 1]
            let coh_norm = clamp01((coh - 0.6) / 1.2);
            // Slow EMA (~5s tau at 44.1kHz/512 block = 86 blocks/s, alpha ≈ 0.002)
            let alpha = clamp(num_samples as f64 / (5.0 * self.sr), 0.0001, 0.05);
            self.coh_slow = (1.0 - alpha) * self.coh_slow + alpha * coh_norm;

            let base_temp = ramped_temp;
            let base_rate = self.params.bilateral_rate;

            if self.coh_slow > 0.7 {
                // High coherence: converge — reduce temperature, slow bilateral
                let factor = (self.coh_slow - 0.7) / 0.3;
                self.effective_temp = clamp(base_temp * (1.0 - 0.15 * factor), 0.10, 0.50);
                self.effective_bilateral_rate = clamp(base_rate * (1.0 - 0.30 * factor), 0.3, 6.0);
            } else if self.coh_slow < 0.3 {
                // Low coherence: explore — increase temperature and rate
                let factor = (0.3 - self.coh_slow) / 0.3;
                self.effective_temp = clamp(base_temp * (1.0 + 0.10 * factor), 0.10, 0.50);
                self.effective_bilateral_rate = clamp(base_rate * (1.0 + 0.30 * factor), 0.3, 6.0);
            } else {
                self.effective_temp = ramped_temp;
            }
            // Update bilateral rate from feedback
            self.bilateral.set_rate(self.effective_bilateral_rate);
        } else {
            self.effective_temp = ramped_temp;
        }

        // --- 0.1 Hz macro-modulation (heart-brain coherence rhythm) ---
        // OU process at ~0.1 Hz modulates width for breathing-like spatial pulsation.
        // Only active when feedback loop is on (therapeutic presets).
        let block_dt = num_samples as f64 / self.sr;
        if self.params.feedback_on {
            self.ou_macro.step(block_dt, 0.0, &mut self.rng);
        }

        // Phi lattice OU drift (slow evolution of ratio preferences)
        if self.params.phi_lattice_on {
            self.phi_lattice.step_drift(block_dt);
        }

        // Weighted-average pan from previous block's grains (for DVF near-field)
        let avg_pan = if self.block_pan_weight > 1e-9 {
            clamp(self.block_pan_sum / self.block_pan_weight, -1.0, 1.0)
        } else {
            0.0
        };
        self.block_pan_sum = 0.0;
        self.block_pan_weight = 0.0;

        // Update tinnitus notch filter (only recalcs on param change)
        self.tinnitus.set_params(self.params.tinnitus_notch_hz, self.params.tinnitus_notch_q, self.sr);

        // Pre-calculate constants
        let lat_inc = clamp(self.params.lat_rate, 1.0, 2000.0) / self.sr;
        let itd_scale = self.params.itd_us * 1.0e-6 * self.sr;
        let ext_cfg = ExternalProcessor::prepare(self.params.externalization, self.sr);

        // T6.3: Coherence-driven spatial morphing — pre-compute modulation
        let mut coh_spatial_mod = if self.params.coherence_spatial {
            // High coherence → wider separation (emphasize bilateral effect)
            let coh = self.dialogue.coherence();
            let coh_norm = clamp01((coh - 0.6) / 1.2);
            // Scale: 0.8 at low coherence, 1.3 at high coherence
            0.8 + 0.5 * coh_norm
        } else {
            1.0
        };

        // Macro-modulation: 0.1 Hz OU scales width ±15% for heart-brain coherence rhythm
        if self.params.feedback_on {
            let macro_norm = 0.5 + 0.5 * self.ou_macro.y.tanh(); // [0, 1]
            coh_spatial_mod *= 0.85 + 0.30 * macro_norm; // [0.85, 1.15]
        }

        // Stochastic Resonance: adapt noise gain from coherence (Collins 1995)
        if self.stoch_res.is_active() {
            let coh = self.dialogue.coherence();
            self.stoch_res.adapt(coh, num_samples, self.sr);
        }

        // Polyrhythm tick for this block
        let _poly_result = if self.params.polyrhythm_on {
            self.polyrhythm.tick(self.sr, num_samples as f64)
        } else {
            (false, false, false)
        };
        let poly_pan_offset = if self.params.polyrhythm_on {
            self.polyrhythm.pan_offset(self.params.polyrhythm_amount)
        } else {
            0.0
        };

        // Modal contralateral spatial mirror — set once per block
        if self.params.modal_on && self.params.modal_contralateral > 1e-6 {
            let bw = if self.burst_engine.enabled {
                self.burst_engine.compute_weight(&self.hawkes)
            } else { 0.0 };
            self.modal_engine.set_contralateral(
                self.burst_centroid,
                bw,
                self.params.modal_contralateral,
            );
        } else {
            self.modal_engine.set_contralateral(0.0, 0.0, 0.0);
        }

        for n in 0..num_samples {
            // Update counters
            self.gap_elapsed += 1;
            
            // Update stochastic processes
            if self.params.thermo || self.params.lattice {
                self.lat_phase += lat_inc;
                if self.lat_phase >= 1.0 {
                    let k = self.lat_phase.floor() as i32;
                    let dt = k as f64 / self.params.lat_rate.max(1.0);
                    
                    if self.params.lattice {
                        for _ in 0..k {
                            self.lattice.step(&mut self.rng);
                        }
                    }
                    
                    if self.params.thermo {
                        let t = clamp01(self.effective_temp);
                        self.ou_pan.sigma = 0.40 * t;
                        self.ou_itd.sigma = 0.35 * t;
                        self.ou_amp.sigma = 0.30 * t;
                        self.ou_rate.sigma = 0.25 * t;
                        self.ou_pan.step(dt, 0.0, &mut self.rng);
                        self.ou_itd.step(dt, 0.0, &mut self.rng);
                        self.ou_amp.step(dt, 0.0, &mut self.rng);
                        self.ou_rate.step(dt, 0.0, &mut self.rng);
                    }
                    
                    if self.params.burst {
                        self.hawkes.tick(dt, &mut self.rng);
                    }
                    
                    self.lat_phase -= k as f64;
                }
            }
            
            // Generate noise: spectral tilt for classic modes (0-2), NoiseGen for extended (3-5)
            let mut nz = if self.params.noise_mode <= 2 {
                // Continuous spectral tilt replaces discrete noise_color
                let white = self.rng.uni_pm1();
                self.spectral_tilt.process(white, self.params.noise_slope)
            } else {
                self.noise_gen.next_sample(&mut self.rng, self.sr)
            };
            // Stochastic resonance gain (Collins 1995: ±2.5 dB around unity)
            if self.stoch_res.is_active() {
                nz *= self.stoch_res.noise_gain();
            }
            // Tinnitus notch filter
            nz = self.tinnitus.process(nz);
            // No pre-clip: output stage soft_tanh handles clipping.
            // Pre-clip was compressing brown noise peaks, flattening spectrum.
            self.ring.write(nz);
            let wi = self.ring.get_write_index();
            
            // Schedule new grain
            self.samples_to_next -= 1;
            if self.samples_to_next <= 0 {
                self.spawn_grain(wi, itd_scale, coh_spatial_mod, poly_pan_offset);
                self.samples_to_next = self.schedule_gap_samples();
            }
            
            // Process all active grains
            let (mut y_l, mut y_r) = (0.0, 0.0);
            
            for grain in self.grains.iter_mut() {
                if !grain.on { continue; }
                if grain.age >= grain.dur {
                    grain.on = false;
                    continue;
                }
                
                let phase = grain.phase();
                let env = Envelope::eval(phase, &grain.env, self.params.envelope_shape);
                
                // Read from ring with ITD + per-grain offset for decorrelation
                let itd = grain.itd;
                let base = (wi + RING_SIZE - grain.ring_offset) & RING_MASK;
                let (mut s_l, mut s_r) = self.ring.read_stereo_itd(base, itd);

                // Per-grain spectral tilt (lightweight single-pole coloring)
                let tilt_diff = grain.grain_tilt - self.params.noise_slope;
                if tilt_diff.abs() > 0.1 {
                    let alpha = 0.997_f64.powf(1.0 + tilt_diff.abs());
                    grain.tilt_z = alpha * grain.tilt_z + (1.0 - alpha) * s_l;
                    let tilt_s = if tilt_diff < 0.0 {
                        grain.tilt_z
                    } else {
                        s_l + 0.3 * (s_l - grain.tilt_z)
                    };
                    let ratio = if s_l.abs() > 1e-12 { tilt_s / s_l } else { 1.0 };
                    s_l = tilt_s;
                    s_r *= ratio;
                }

                // Apply pan and ILD
                s_l *= grain.pan_l * grain.gain_l;
                s_r *= grain.pan_r * grain.gain_r;
                
                // Crossfeed
                if grain.crossfeed.abs() > 1.0e-6 {
                    let base_l = s_l;
                    let base_r = s_r;
                    s_l = base_l + grain.crossfeed * base_r;
                    s_r = base_r + grain.crossfeed * base_l;
                }
                
                // Sample rate crush (sample and hold)
                grain.sr_hold_cnt -= 1;
                if grain.sr_hold_cnt <= 0 {
                    grain.held_l = s_l;
                    grain.held_r = s_r;
                    grain.sr_hold_cnt = grain.sr_hold_n;
                }
                s_l = grain.held_l;
                s_r = grain.held_r;
                
                // Bit crush
                if grain.q_levels > 0 {
                    let q = grain.q_levels as f64;
                    s_l = (s_l * q).round() / q;
                    s_r = (s_r * q).round() / q;
                }
                
                // Glitch effects
                match grain.kind {
                    GrainKind::Stutter => {
                        if (grain.age & 7) == 0 {
                            s_l *= 0.2;
                            s_r *= 0.2;
                        }
                    }
                    _ => {}
                }
                
                // IPD (allpass decorrelation)
                if grain.ipd_coeff.abs() > 1.0e-6 {
                    let (new_l, new_z_l) = allpass(s_l, grain.ipd_coeff, grain.ipd_z_l);
                    let (new_r, new_z_r) = allpass(s_r, -grain.ipd_coeff, grain.ipd_z_r);
                    s_l = new_l;
                    s_r = new_r;
                    grain.ipd_z_l = new_z_l;
                    grain.ipd_z_r = new_z_r;
                }
                
                // Head shadow (contralateral LP filter)
                if grain.shadow_a > 1.0e-6 {
                    if grain.shadow_left {
                        let (new_l, new_z) = shadow_lp(s_l, grain.shadow_a, grain.shadow_z_l);
                        s_l = new_l;
                        grain.shadow_z_l = new_z;
                    }
                    if grain.shadow_right {
                        let (new_r, new_z) = shadow_lp(s_r, grain.shadow_a, grain.shadow_z_r);
                        s_r = new_r;
                        grain.shadow_z_r = new_z;
                    }
                }

                // Distance LP (HF rolloff with distance — phi spatial pipeline)
                if grain.dist_lp_a > 1e-6 {
                    grain.dist_lp_z_l = grain.dist_lp_a * grain.dist_lp_z_l
                        + (1.0 - grain.dist_lp_a) * s_l;
                    s_l = grain.dist_lp_z_l;
                    grain.dist_lp_z_r = grain.dist_lp_a * grain.dist_lp_z_r
                        + (1.0 - grain.dist_lp_a) * s_r;
                    s_r = grain.dist_lp_z_r;
                }

                // Accumulate with envelope
                y_l += grain.amp * env * s_l;
                y_r += grain.amp * env * s_r;
                
                grain.age += 1;
            }
            
            // Modal resonator (post grain-sum, pre external)
            if self.params.modal_on {
                let (ml, mr, contra_mono, contra_amt) = self.modal_engine.process(y_l, y_r);
                y_l = ml;
                y_r = mr;

                // Spatialize contralateral signal through phi head model
                if contra_amt > 1e-6 {
                    let mirror_pan = -self.burst_centroid;

                    // 1. Equal-power pan
                    let (mp_l, mp_r) = pan_equal_power(mirror_pan, self.params.width);

                    // 2. ITD from phi head geometry
                    let head = compute_head_result(
                        &self.phi_geom, self.sr, mirror_pan, 0.0, self.params.phi_distance,
                    );
                    // head.itd_samples is already in samples; scale by user ITD param
                    let itd_samples = head.itd_samples
                        * (self.params.itd_us / 700.0);
                    let half_itd = itd_samples * 0.5;

                    // Write contra_mono to delay ring
                    self.contra_delay[self.contra_delay_pos] = contra_mono;
                    self.contra_delay_pos = (self.contra_delay_pos + 1) % 64;

                    // Read with ITD: base delay allows both-direction offsets.
                    // mirror_pan > 0 → source on right → right ear closer → less delay
                    let base_delay = 16.0_f64;
                    let read_l = contra_read_lerp(
                        &self.contra_delay, self.contra_delay_pos,
                        base_delay + half_itd,
                    );
                    let read_r = contra_read_lerp(
                        &self.contra_delay, self.contra_delay_pos,
                        base_delay - half_itd,
                    );

                    // 3. ILD (level difference)
                    let ild_db_val = map_ild_db(self.params.ild_db, mirror_pan, 0.5);
                    let ild_lin = db_to_lin(ild_db_val);
                    // Contralateral ear is attenuated
                    let (gain_l, gain_r) = if mirror_pan > 0.0 {
                        (ild_lin, 1.0)  // source right → left ear attenuated
                    } else {
                        (1.0, ild_lin)  // source left → right ear attenuated
                    };

                    let mut cs_l = read_l * mp_l * gain_l;
                    let mut cs_r = read_r * mp_r * gain_r;

                    // 4. Head shadow (LP on contralateral ear)
                    let shadow_cutoff = 800.0 + 7200.0 * (1.0 - mirror_pan.abs());
                    let shadow_a = (-TWO_PI * shadow_cutoff / self.sr).exp();
                    if mirror_pan > 0.0 {
                        // Source on right → shadow on left
                        let (new_l, new_z) = shadow_lp(cs_l, shadow_a, self.contra_shadow_z_l);
                        cs_l = new_l;
                        self.contra_shadow_z_l = new_z;
                    } else {
                        // Source on left → shadow on right
                        let (new_r, new_z) = shadow_lp(cs_r, shadow_a, self.contra_shadow_z_r);
                        cs_r = new_r;
                        self.contra_shadow_z_r = new_z;
                    }

                    // 5. Crossfeed (subtle inter-ear bleed)
                    let focus = 0.25 + 0.75
                        * (mirror_pan.abs() * std::f64::consts::FRAC_PI_2).sin().powf(1.35);
                    let xfeed = (1.0 - focus) * 0.18;
                    let bl = cs_l;
                    let br = cs_r;
                    cs_l = bl + xfeed * br;
                    cs_r = br + xfeed * bl;

                    // Add spatialized contralateral to output
                    y_l += cs_l * contra_amt;
                    y_r += cs_r * contra_amt;
                }
            }

            // Binaural beat (separate from noise path — Wahbeh 2007)
            if self.params.binaural_on {
                let (bl, br) = self.binaural.process_sample(
                    self.sr,
                    self.params.binaural_carrier_hz,
                    self.params.binaural_beat_hz,
                    self.params.binaural_level,
                );
                y_l += bl;
                y_r += br;
            }

            // Isochronic tone (mono, added to both channels)
            if self.params.isochronic_on {
                let iso = self.isochronic.process_sample(
                    self.sr,
                    self.params.isochronic_carrier_hz,
                    self.params.isochronic_rate_hz,
                    self.params.isochronic_duty,
                    self.params.isochronic_level,
                );
                y_l += iso;
                y_r += iso;
            }

            // DVF near-field (T6.1) — per-ear high-shelf boost at close distance
            if self.dvf.is_active() {
                let (dl, dr) = self.dvf.process(y_l, y_r, avg_pan);
                y_l = dl;
                y_r = dr;
            }

            // External externalisation (block-level cross-channel feedback delay)
            self.external_proc.process_sample(&ext_cfg, &mut y_l, &mut y_r);

            // Room reverb (T6.2) — phi-ratio Schroeder allpass
            let (rl, rr) = self.room.process(y_l, y_r);
            y_l = rl;
            y_r = rr;

            // Body resonance safety notch filters (P0 — chest cavity 5-8 Hz, eyeball 19 Hz)
            y_l = biquad_tick(y_l, &self.safety_notch_chest_coeff, &mut self.safety_notch_chest_zl);
            y_r = biquad_tick(y_r, &self.safety_notch_chest_coeff, &mut self.safety_notch_chest_zr);
            y_l = biquad_tick(y_l, &self.safety_notch_eye_coeff, &mut self.safety_notch_eye_zl);
            y_r = biquad_tick(y_r, &self.safety_notch_eye_coeff, &mut self.safety_notch_eye_zr);

            // Soft clip output
            out_l[n] = soft_tanh(y_l * OUT_DRIVE) / OUT_DRIVE;
            out_r[n] = soft_tanh(y_r * OUT_DRIVE) / OUT_DRIVE;
            
            self.sample_counter += 1;
        }
    }
    
    fn spawn_grain(&mut self, wi: usize, itd_scale: f64, coh_spatial_mod: f64, poly_pan_offset: f64) {
        let gi = self.grains.find_free();
        if gi < 0 { return; }
        
        let grain = match self.grains.get_mut(gi as usize) {
            Some(g) => g,
            None => return,
        };
        grain.reset();
        
        // Generate quasi-random values
        let u1 = self.w_phi.next();
        let u2 = self.w_s2.next();
        let u3 = self.w_pl.next();
        let u4 = self.rng.uni01();
        let u5 = self.rng.uni01();
        let u6 = self.rng.uni01();
        
        // Hemisphere coupling
        let gap_samples = self.gap_elapsed as f64;
        let prev_dur = self.last_dur_samples.max(1.0);
        let ratio_gap = gap_samples / prev_dur;
        let coupling = clamp01(self.params.hemis_coupling);
        // φ-corrected coefficients
        let time_weight = clamp(
            HEMI_WEIGHT_MIN + HEMI_WEIGHT_MAX * (ratio_gap / (ratio_gap + 1.0)),
            HEMI_WEIGHT_MIN,
            1.0
        );
        let hemi = coupling * time_weight;
        
        // Lattice modulation
        let lat_u = if self.params.lattice {
            let v = self.lattice.probe(self.w_phi.next());
            self.lat_last_v = v;
            0.5 + 0.5 * v.tanh()
        } else {
            0.5
        };
        
        // OU modulation
        let oup = if self.params.thermo { clamp(self.ou_pan.y, -1.0, 1.0) } else { 0.0 };
        let oua = if self.params.thermo { 
            map_phi_range(INV_PHI, PHI, 0.5 + 0.5 * self.ou_amp.y.tanh()) 
        } else { 
            1.0 
        };
        let oui = if self.params.thermo { self.ou_itd.y } else { 0.0 };
        
        // Calculate grain parameters
        let amp_shape = u2.max(1e-9).powf(0.35);
        let amp_lat = map_phi_range(INV_PHI, PHI, lat_u);
        grain.amp = AMP_NORM * amp_shape * amp_lat * oua;
        
        let mut pan = 2.0 * u3 - 1.0;
        if self.params.lattice {
            pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
        }
        pan = clamp((1.0 - hemi) * pan - hemi * self.prev_pan, -1.0, 1.0);
        grain.pan = pan;
        
        // Duration
        let base = clamp(self.params.baselen_ms, MIN_BASE_LENGTH_MS, 2000.0) * 0.001 * self.sr;
        let kexp = (2.0 * u1 - 1.0) * clamp01(self.params.len_phi);
        let len = clamp(base * PHI.powf(kexp), MIN_GRAIN_SAMPLES, self.sr * 4.0);
        grain.dur = len as u32;

        // --- Dialogue evaluation FIRST (interhemispheric coherence) ---
        // Dialogue must see the raw stochastic pan to detect phi-ratio handshakes.
        // Previously burst modified pan before dialogue, destroying phi relationships.
        let burst_weight = if self.burst_engine.enabled {
            self.burst_engine.compute_weight(&self.hawkes)
        } else {
            0.0
        };
        let dialogue_params = DialogueParams {
            strength: self.params.dialogue_strength,
            memory: self.params.dialogue_memory,
            phi_mix: self.params.dialogue_phi_mix,
            enabled: self.params.dialogue_on,
            burst_weight,
        };
        let dial_result = self.dialogue.evaluate(
            &dialogue_params, pan, grain.amp, len, gap_samples,
        );
        if self.params.dialogue_on {
            pan = dial_result.pan;
            grain.pan = pan;
            grain.amp *= dial_result.amp_scale;
            let new_len = clamp(len * dial_result.dur_scale, MIN_GRAIN_SAMPLES, self.sr * 4.0);
            grain.dur = new_len as u32;
        }
        self.dialogue.commit(&dialogue_params, &dial_result, true);

        // --- Burst position modulation AFTER dialogue ---
        // Burst's center_pull is spatial clustering, independent of phi-ratio detection.
        if self.burst_engine.enabled {
            let dur_norm = clamp01(len / (base * PHI));
            let expected_gap = if self.params.rate > 1.0e-6 {
                self.sr / self.params.rate
            } else {
                self.sr * 0.25
            };
            let gap_norm = clamp01(gap_samples / expected_gap);
            let br = self.burst_engine.apply_position(&self.hawkes, pan, dur_norm, gap_norm);
            pan = br.pan;
            grain.pan = pan;
            grain.amp *= br.amp_scale;
        }

        // Update burst centroid EMA (tracks cluster center of mass)
        if self.burst_engine.enabled {
            let bw = self.burst_engine.compute_weight(&self.hawkes);
            if bw > 0.01 {
                let alpha = clamp(self.burst_centroid_alpha * bw, 0.0001, 0.1);
                self.burst_centroid = (1.0 - alpha) * self.burst_centroid + alpha * pan;
            } else {
                // Decay toward center when no burst
                self.burst_centroid *= 0.999;
            }
        }

        // --- Phi-Pan (phi-ratio alternation) ---
        if self.params.phi_pan {
            pan = self.phi_pan_proc.next(pan);
            grain.pan = pan;
        }

        // --- Bilateral oscillator (EMDR-style deterministic L-R sweep) ---
        // T0.1 fix: pass gap_samples so phase advances by actual elapsed time
        if self.params.bilateral_on {
            pan = self.bilateral.apply(pan, self.sr, gap_samples);
            grain.pan = pan;
        }

        // --- Polyrhythm pan offset (T5.5) ---
        // p_pulse pushes left, q_pulse pushes right, coincidence → center
        if self.params.polyrhythm_on && poly_pan_offset.abs() > 1e-6 {
            pan = clamp(pan + poly_pan_offset, -1.0, 1.0);
            grain.pan = pan;
        }

        // --- Inter-grain decorrelation: quasi-random ring offset ---
        // Each grain reads from a different region of the ring buffer.
        // Offset capped at 4096 samples (~93ms at 44.1kHz) for fresh content.
        grain.ring_offset = (self.w_s2.next() * 4096.0) as usize;

        // Accumulate weighted pan for DVF near-field (block-level average)
        self.block_pan_sum += pan * grain.amp;
        self.block_pan_weight += grain.amp;

        // Binaural — T6.3: coherence-driven spatial morphing scales width
        let eff_width = self.params.width * coh_spatial_mod;
        let (pan_l, pan_r) = pan_equal_power(pan, eff_width);
        grain.pan_l = pan_l;
        grain.pan_r = pan_r;
        
        // ITD: phi-model ellipsoid Woodworth (direction-dependent head radius)
        let max_itd_samples = clamp(self.params.itd_us, 0.0, 1600.0) * 1.0e-6 * self.sr;
        let mut itd = if max_itd_samples > 1e-9 {
            let phi_head = compute_head_result(&self.phi_geom, self.sr, pan, 0.0, self.params.phi_distance);
            let scale = max_itd_samples / self.phi_itd_max;
            let base_itd = phi_head.itd_samples * scale;
            // Stochastic jitter (±8% of max)
            let jitter = (2.0 * clamp01(u2) - 1.0) * 0.08 * max_itd_samples;
            clamp(base_itd + jitter, -max_itd_samples, max_itd_samples)
        } else {
            0.0
        };
        if self.params.lattice {
            itd += ((2.0 * lat_u - 1.0) * 0.33 + 0.33 * oui) * itd_scale;
        }
        grain.itd = (1.0 - hemi) * itd - hemi * self.prev_itd;
        
        let ild_db = map_ild_db(self.params.ild_db, pan, u5);
        grain.gain_l = db_to_lin(ild_db);
        grain.gain_r = db_to_lin(-ild_db);
        
        // Focus and crossfeed
        let lateral = (pan * PI * 0.5).sin();
        grain.focus = 0.25 + 0.75 * lateral.abs().powf(1.35);
        grain.crossfeed = (1.0 - grain.focus) * 0.18;
        
        // IPD
        let ipd_amt = clamp01(self.params.spat_ipd);
        if ipd_amt > 1.0e-6 {
            let ipd_shape = (2.0 * u6 - 1.0).abs();
            let ipd_base = 0.18 + 0.55 * ipd_amt;
            let ipd_spread = 0.25 * ipd_amt;
            grain.ipd_coeff = clamp(ipd_base + ipd_spread * (ipd_shape - 0.5), 0.0, 0.95);
            grain.ipd_coeff *= 0.7 + 0.3 * grain.focus;
        }
        
        // Head shadow
        let shadow_amt = clamp01(self.params.spat_shadow);
        let shadow_factor = shadow_amt * pan.abs();
        if shadow_factor > 1.0e-6 {
            let fc_min: f64 = 800.0;
            let fc_max: f64 = fc_min.max((0.45 * self.sr).min(8000.0));
            let fc = map_phi_range(fc_min, fc_max, 1.0 - shadow_factor);
            grain.shadow_a = clamp((-TWO_PI * fc / self.sr).exp(), 0.0, 0.9999);
            grain.shadow_left = pan > 0.0;
            grain.shadow_right = pan < 0.0;
        }

        // --- Phi spatial pipeline: pinna + distance + air absorption ---
        // Pinna amplitude modulation (azimuth-dependent gain from pinna geometry).
        // At pan=0 (center), amp_scale≈0.35; at pan=1.0 (lateral), amp_scale≈1.0.
        // This creates a ~16% dip at center during bilateral sweep (at spat_pinna=0.25),
        // which reinforces the bilateral percept. At higher values the dip deepens.
        let pinna_amt = clamp01(self.params.spat_pinna);
        if pinna_amt > 1e-6 {
            let pinna_tuning = PinnaTuning {
                elev_norm: clamp01(self.params.phi_elev),
                ..PinnaTuning::default()
            };
            let pinna_resp = design_pinna_response(&self.phi_geom, pan, &pinna_tuning);
            // Mix: at 0 no effect, at 1 full pinna attenuation
            grain.amp *= 1.0 - pinna_amt * (1.0 - pinna_resp.amp_scale);
        }

        // Distance response: gain attenuation + LP HF rolloff + air absorption.
        // phi_distance is clamped to [0,1] here because compute_distance_response maps
        // to physical distance [0.8m, 4.0m]. The ITD path above uses unclamped phi_distance
        // (default 1.5) where values >1 reduce ITD (smaller angle subtended at distance).
        let dist_amt = clamp01(self.params.spat_distance);
        if dist_amt > 1e-6 {
            let phi_dist = clamp01(self.params.phi_distance);
            let dist_resp = compute_distance_response(self.sr, phi_dist);
            // Mix: at 0 no attenuation, at 1 full distance model
            grain.amp *= 1.0 - dist_amt * (1.0 - dist_resp.direct_gain);
            // LP coefficient scaled by distance amount
            grain.dist_lp_a = dist_resp.lowpass_alpha * dist_amt;

            // Air absorption at reference 4kHz (ISO 9613-1)
            let distance_m = 0.8 + 3.2 * phi_dist;
            let air_atten = compute_air_absorption(4000.0, distance_m);
            grain.amp *= 1.0 - dist_amt * (1.0 - air_atten);
        }

        // Glitch kind
        grain.kind = GrainKind::choose(self.params.glitch_mix, u4);
        
        // SR/bit crush
        grain.sr_hold_n = map_sr_hold_base(self.params.srcrush_amt, u1);
        grain.sr_hold_cnt = grain.sr_hold_n;
        grain.q_levels = (1 << (map_bits(self.params.bitcrush_amt) - 1)) - 1;
        
        // Envelope — phi lattice personality or global params
        let (env_a, env_d, env_s, env_r) = if self.params.phi_lattice_on && self.params.phi_personality > 1e-6 {
            let (la, ld, ls, lr) = self.phi_lattice.grain_adsr();
            let p = clamp01(self.params.phi_personality);
            (
                (1.0 - p) * self.params.env_attack + p * la,
                (1.0 - p) * self.params.env_decay + p * ld,
                (1.0 - p) * self.params.env_sustain + p * ls,
                (1.0 - p) * self.params.env_release + p * lr,
            )
        } else {
            (self.params.env_attack, self.params.env_decay, self.params.env_sustain, self.params.env_release)
        };

        // Store personality in grain
        grain.grain_attack = env_a;
        grain.grain_decay = env_d;
        grain.grain_sustain = env_s;
        grain.grain_release = env_r;

        // Bilateral safety: clamp attack so onset < 5ms for CC stimulation
        let mut final_env_a = env_a;
        let dur_ms = len * 1000.0 / self.sr;
        if self.params.bilateral_on && final_env_a * dur_ms > 5.0 {
            final_env_a = 5.0 / dur_ms;
            grain.grain_attack = final_env_a;
        }

        grain.env = self.envelope.make_shape(
            final_env_a, env_d, env_s, env_r,
            gap_samples, len, pan.abs(),
        );

        // Per-grain spectral tilt from phi lattice
        if self.params.phi_lattice_on && self.params.phi_personality > 1e-6 {
            grain.grain_tilt = self.phi_lattice.grain_spectral_tilt(self.params.noise_slope);
        } else {
            grain.grain_tilt = self.params.noise_slope;
        }

        // Activate grain
        grain.on = true;
        grain.age = 0;
        
        // Update previous state
        self.prev_pan = pan;
        self.prev_itd = grain.itd;
        self.prev_ild = ild_db;
        self.last_gap_samples = gap_samples;
        self.last_dur_samples = len;
        self.gap_elapsed = 0;
    }
    
    fn schedule_gap_samples(&mut self) -> i32 {
        // T1.2: After handshake, use pre-computed Fibonacci-ratio gaps ("temporal echo")
        if let Some(queued) = self.dialogue.pop_queued_gap() {
            return queued.round().max(1.0) as i32;
        }

        let mut rate = clamp(self.params.rate, 0.0, MAX_EVENT_RATE_HZ);

        // Thermo OU rate modulation — apply BEFORE nesting so nesting quantizes the final rate
        if self.params.thermo {
            let ur = 0.5 + 0.5 * self.ou_rate.y.tanh();
            let rate_phi = map_phi_range(
                (self.params.rate / PHI).max(0.001),
                self.params.rate * PHI,
                ur
            );
            rate = clamp(rate_phi, 0.0, MAX_EVENT_RATE_HZ);
        }

        // T2.3: Theta-gamma nesting — quantize grain rate to integer multiple of bilateral rate
        // Lisman-Jensen 2013: 4-8 gamma cycles per theta cycle.
        // Must run AFTER thermo so the OU-modulated rate gets quantized, not overwritten.
        if self.params.bilateral_nesting && self.params.bilateral_on && self.params.bilateral_rate > 0.1 {
            let bi_rate = self.params.bilateral_rate;
            let ratio = (rate / bi_rate).round().max(4.0).min(8.0); // 4:1 to 8:1
            rate = bi_rate * ratio;
        }
        
        if rate <= 1e-6 {
            return (self.sr * 0.25).max(1.0) as i32;
        }
        
        let t = self.sample_counter as f64 / self.sr;
        let mut lambda = rate * (1.0 + 0.2 * (TWO_PI * (t * INV_PHI)).sin());
        lambda = lambda.max(1e-3);
        
        if self.params.burst {
            lambda += 0.3 * self.hawkes.lambda;
        }
        
        let u = self.rng.uni01().max(1.0e-12);
        let mut gap_sec = -u.ln() / lambda;
        
        let base_rate = clamp(self.params.rate, 0.0, MAX_EVENT_RATE_HZ);
        let cap_sec = if base_rate > 1e-6 {
            30.0_f64.min(4.0 / base_rate)
        } else {
            0.25
        };
        gap_sec = gap_sec.min(cap_sec);
        
        (gap_sec * self.sr).round().max(1.0) as i32
    }
}

/// Allpass filter for IPD decorrelation
#[inline]
fn allpass(x: f64, a: f64, z: f64) -> (f64, f64) {
    let y = -a * x + z;
    let new_z = x + a * y;
    (y, new_z)
}

/// Simple lowpass for head shadow
#[inline]
fn shadow_lp(x: f64, a: f64, z: f64) -> (f64, f64) {
    let y = (1.0 - a) * x + a * z;
    (y, y)
}

/// Compute 2nd-order biquad notch (band-reject) filter coefficients.
/// Returns [b0, b1, b2, a1, a2].
/// `center_hz`: notch center frequency
/// `sr`: sample rate
/// `q`: quality factor (bandwidth control, higher = narrower notch)
/// Attenuation at center frequency is theoretically infinite for ideal notch;
/// practical Q of 1.0–2.0 gives 12–20 dB rejection over the target band.
#[inline]
fn compute_notch_coeffs(center_hz: f64, sr: f64, q: f64) -> [f64; 5] {
    let w0 = TWO_PI * center_hz / sr;
    let alpha = w0.sin() / (2.0 * q);
    let cos_w0 = w0.cos();

    let b0 = 1.0;
    let b1 = -2.0 * cos_w0;
    let b2 = 1.0;
    let a0 = 1.0 + alpha;
    let a1 = -2.0 * cos_w0;
    let a2 = 1.0 - alpha;

    // Normalize by a0
    [b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0]
}

/// Apply biquad filter (transposed direct form II) to a single sample.
/// `c`: [b0, b1, b2, a1, a2], `z`: [z1, z2] (state, mutated in place).
/// Returns filtered sample.
#[inline]
fn biquad_tick(x: f64, c: &[f64; 5], z: &mut [f64; 2]) -> f64 {
    let y = c[0] * x + z[0];
    z[0] = c[1] * x - c[3] * y + z[1];
    z[1] = c[2] * x - c[4] * y;
    y
}

/// Read from the 64-sample contralateral delay ring with linear interpolation.
/// `write_pos` is the next-write position (most recent sample is at write_pos - 1).
/// `delay` is in fractional samples (clamped to [0, 62]).
#[inline]
fn contra_read_lerp(buf: &[f64; 64], write_pos: usize, delay: f64) -> f64 {
    let d = delay.clamp(0.0, 62.0);
    let i = d as usize;
    let frac = d - i as f64;
    let idx0 = (write_pos + 64 - 1 - i) % 64;
    let idx1 = (write_pos + 64 - 2 - i) % 64;
    buf[idx0] * (1.0 - frac) + buf[idx1] * frac
}

/// Python module initialization
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Rng>()?;
    m.add_class::<Weyl>()?;
    m.add_class::<NoiseColor>()?;
    m.add_class::<NoiseColorState>()?;
    m.add_class::<NoiseMode>()?;
    m.add_class::<NoiseGen>()?;
    m.add_class::<GrainKind>()?;
    m.add_class::<RingBuffer>()?;
    m.add_class::<OrnsteinUhlenbeck>()?;
    m.add_class::<Lattice>()?;
    m.add_class::<Hawkes>()?;
    m.add_class::<EnvelopeShape>()?;
    m.add_class::<Envelope>()?;
    m.add_class::<Grain>()?;
    m.add_class::<GrainPool>()?;
    m.add_class::<Params>()?;
    m.add_class::<Engine>()?;
    m.add_class::<PhiModel>()?;

    // Export constants
    m.add("PHI", PHI)?;
    m.add("INV_PHI", INV_PHI)?;
    m.add("PHI_SQ", PHI_SQ)?;
    m.add("INV_PHI_SQ", INV_PHI_SQ)?;
    m.add("INV_PHI_CU", INV_PHI_CU)?;
    m.add("MAX_GRAINS", MAX_GRAINS)?;
    m.add("RING_SIZE", RING_SIZE)?;

    Ok(())
}
