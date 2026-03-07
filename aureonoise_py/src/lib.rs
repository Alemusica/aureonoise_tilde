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
mod external;
mod modal;
mod dialogue;
mod binaural;
mod isochronic;
mod tinnitus;
mod dvf;
mod room;
mod polyrhythm;

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
pub use phi_model::{PhiModel, GeometryConfig, Geometry, build_geometry, compute_head_result};
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
    
    // Timbre
    #[pyo3(get, set)]
    pub noise_color: i32,
    #[pyo3(get, set)]
    pub color_amt: f64,
    #[pyo3(get, set)]
    pub vhs_wow: f64,
    #[pyo3(get, set)]
    pub vhs_flutter: f64,
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
            
            // Timbre
            noise_color: 1, // Pink
            color_amt: 0.65,
            vhs_wow: 0.35,
            vhs_flutter: 0.25,
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

    // Coherence feedback loop (BAC-inspired)
    coh_slow: f64,           // slow EMA of coherence (~5s tau)
    effective_temp: f64,     // feedback-modulated temperature
    effective_bilateral_rate: f64, // feedback-modulated bilateral rate

    // Temperature ramp
    temp_ramp_samples: f64,  // total ramp duration in samples (0 = instant)
    temp_ramp_elapsed: f64,  // samples elapsed in current ramp
    temp_start: f64,         // starting temperature
    temp_target: f64,        // target temperature

    // LFO
    lfo_wow_phase: f64,
    lfo_flut_phase: f64,

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
            ou_pan: OrnsteinUhlenbeck::new(0.60, 0.0),
            ou_itd: OrnsteinUhlenbeck::new(0.40, 0.0),
            ou_amp: OrnsteinUhlenbeck::new(0.80, 0.0),
            ou_rate: OrnsteinUhlenbeck::new(1.20, 0.0),
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
            lfo_wow_phase: 0.0,
            lfo_flut_phase: 0.0,
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
        self.external_proc.reset();
        self.dialogue.reset();
        self.phi_pan_proc.reset();
        self.bilateral.reset();
        self.modal_engine.reset();
        self.dvf.reset();
        self.room.reset();
        self.polyrhythm.reset();
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
        self.lfo_wow_phase = 0.0;
        self.lfo_flut_phase = 0.0;
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
    }
}

impl Engine {
    /// Internal block processing
    pub fn process_block(&mut self, out_l: &mut [f64], out_r: &mut [f64]) {
        let num_samples = out_l.len().min(out_r.len());

        // --- Temperature ramp (T4.2) ---
        if self.temp_ramp_samples > 0.0 && self.temp_ramp_elapsed < self.temp_ramp_samples {
            self.temp_ramp_elapsed += num_samples as f64;
            let t = clamp01(self.temp_ramp_elapsed / self.temp_ramp_samples);
            self.effective_temp = self.temp_start + t * (self.temp_target - self.temp_start);
        }

        // --- Coherence feedback loop (T4.1, BAC-inspired) ---
        if self.params.feedback_on {
            let coh = self.dialogue.coherence();
            // Normalize: coherence typically in [0.6, 1.8], map to [0, 1]
            let coh_norm = clamp01((coh - 0.6) / 1.2);
            // Slow EMA (~5s tau at 44.1kHz/512 block = 86 blocks/s, alpha ≈ 0.002)
            let alpha = clamp(num_samples as f64 / (5.0 * self.sr), 0.0001, 0.05);
            self.coh_slow = (1.0 - alpha) * self.coh_slow + alpha * coh_norm;

            let base_temp = self.effective_temp;
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
            }
            // Update bilateral rate from feedback
            self.bilateral.set_rate(self.effective_bilateral_rate);
        }

        // Update tinnitus notch filter (only recalcs on param change)
        self.tinnitus.set_params(self.params.tinnitus_notch_hz, self.params.tinnitus_notch_q, self.sr);

        // Pre-calculate constants
        let wow_hz = map_phi_range(0.1, 1.5, clamp01(self.params.vhs_wow));
        let flt_hz = map_phi_range(7.0, 12.0, clamp01(self.params.vhs_flutter));
        let inc_wow = wow_hz / self.sr;
        let inc_flt = flt_hz / self.sr;
        let lat_inc = clamp(self.params.lat_rate, 1.0, 2000.0) / self.sr;
        let itd_scale = self.params.itd_us * 1.0e-6 * self.sr;
        let ext_cfg = ExternalProcessor::prepare(self.params.externalization, self.sr);

        // T6.3: Coherence-driven spatial morphing — pre-compute modulation
        let coh_spatial_mod = if self.params.coherence_spatial {
            // High coherence → wider separation (emphasize bilateral effect)
            let coh = self.dialogue.coherence();
            let coh_norm = clamp01((coh - 0.6) / 1.2);
            // Scale: 0.8 at low coherence, 1.3 at high coherence
            0.8 + 0.5 * coh_norm
        } else {
            1.0
        };

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

        for n in 0..num_samples {
            // Update counters
            self.gap_elapsed += 1;
            
            // Update LFOs
            self.lfo_wow_phase += inc_wow;
            if self.lfo_wow_phase >= 1.0 { self.lfo_wow_phase -= 1.0; }
            self.lfo_flut_phase += inc_flt;
            if self.lfo_flut_phase >= 1.0 { self.lfo_flut_phase -= 1.0; }
            
            let wow = (TWO_PI * self.lfo_wow_phase).sin();
            let flt = (TWO_PI * self.lfo_flut_phase).sin();
            let vhs_mod = 0.5 * wow + 0.5 * flt;
            
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
            
            // Generate and write noise to ring buffer (dispatches to all 6 modes)
            let mut nz = self.noise_gen.next_sample(&mut self.rng, self.sr);
            // Continuous spectral tilt (overrides discrete noise_color for classic modes)
            if self.params.noise_mode <= 2 {
                let white = self.rng.uni_pm1();
                nz = self.spectral_tilt.process(white, self.params.noise_slope);
            }
            // Tinnitus notch filter
            nz = self.tinnitus.process(nz);
            nz = soft_tanh(nz * 1.2);
            self.ring.write(nz);
            let wi = self.ring.get_write_index();
            
            // Schedule new grain
            self.samples_to_next -= 1;
            if self.samples_to_next <= 0 {
                self.spawn_grain(wi, vhs_mod, itd_scale, coh_spatial_mod, poly_pan_offset);
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
                let env = Envelope::eval(phase, &grain.env);
                
                // Read from ring with ITD + per-grain offset for decorrelation
                let itd = grain.itd + vhs_mod * 0.25 * itd_scale;
                let base = (wi + RING_SIZE - grain.ring_offset) & RING_MASK;
                let (mut s_l, mut s_r) = self.ring.read_stereo_itd(base, itd);
                
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
                    GrainKind::VhsDrop => {
                        let att = 0.5 + 0.5 * (1.0 - vhs_mod.abs());
                        s_l *= att;
                        s_r *= att;
                    }
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
                
                // Accumulate with envelope
                y_l += grain.amp * env * s_l;
                y_r += grain.amp * env * s_r;
                
                grain.age += 1;
            }
            
            // Modal resonator (post grain-sum, pre external)
            if self.params.modal_on {
                let (ml, mr) = self.modal_engine.process(y_l, y_r);
                y_l = ml;
                y_r = mr;
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
                let (dl, dr) = self.dvf.process(y_l, y_r, 0.0);
                y_l = dl;
                y_r = dr;
            }

            // External externalisation (block-level cross-channel feedback delay)
            self.external_proc.process_sample(&ext_cfg, &mut y_l, &mut y_r);

            // Room reverb (T6.2) — phi-ratio Schroeder allpass
            let (rl, rr) = self.room.process(y_l, y_r);
            y_l = rl;
            y_r = rr;

            // Soft clip output
            out_l[n] = soft_tanh(y_l * OUT_DRIVE) / OUT_DRIVE;
            out_r[n] = soft_tanh(y_r * OUT_DRIVE) / OUT_DRIVE;
            
            self.sample_counter += 1;
        }
    }
    
    fn spawn_grain(&mut self, wi: usize, vhs_mod: f64, itd_scale: f64, coh_spatial_mod: f64, poly_pan_offset: f64) {
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

        // Binaural — T6.3: coherence-driven spatial morphing scales width
        let eff_width = self.params.width * coh_spatial_mod;
        let (pan_l, pan_r) = pan_equal_power(pan, eff_width);
        grain.pan_l = pan_l;
        grain.pan_r = pan_r;
        
        // ITD: phi-model ellipsoid Woodworth (direction-dependent head radius)
        let max_itd_samples = clamp(self.params.itd_us, 0.0, 1600.0) * 1.0e-6 * self.sr;
        let mut itd = if max_itd_samples > 1e-9 {
            let phi_head = compute_head_result(&self.phi_geom, self.sr, pan, 0.0, 1.0);
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
        
        // Glitch kind
        grain.kind = GrainKind::choose(self.params.glitch_mix, u4);
        
        // SR/bit crush
        grain.sr_hold_n = map_sr_hold_base(self.params.srcrush_amt, u1);
        grain.sr_hold_cnt = grain.sr_hold_n;
        grain.q_levels = (1 << (map_bits(self.params.bitcrush_amt) - 1)) - 1;
        
        // Envelope
        grain.env = self.envelope.make_shape(
            self.params.env_attack,
            self.params.env_decay,
            self.params.env_sustain,
            self.params.env_release,
            gap_samples,
            len,
            pan.abs(),
        );
        
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

        // T2.3: Theta-gamma nesting — quantize grain rate to integer multiple of bilateral rate
        if self.params.bilateral_nesting && self.params.bilateral_on && self.params.bilateral_rate > 0.1 {
            let bi_rate = self.params.bilateral_rate;
            let ratio = (rate / bi_rate).round().max(4.0).min(8.0); // 4:1 to 8:1
            rate = bi_rate * ratio;
        }

        if self.params.thermo {
            let ur = 0.5 + 0.5 * self.ou_rate.y.tanh();
            let rate_phi = map_phi_range(
                (self.params.rate / PHI).max(0.001),
                self.params.rate * PHI,
                ur
            );
            rate = clamp(rate_phi, 0.0, MAX_EVENT_RATE_HZ);
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
