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

pub use constants::*;
pub use rng::Rng;
pub use weyl::Weyl;
pub use math::*;
pub use noise::{NoiseColor, NoiseColorState, GrainKind};
pub use ring::RingBuffer;
pub use stoch::{OrnsteinUhlenbeck, Lattice, Hawkes};
pub use envelope::{Envelope, EnvelopeShape};
pub use grain::{Grain, GrainPool};

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
            temperature: 0.45,
            lat_rate: 250.0,
            lat_eps: INV_PHI_CU,
            lat_gamma: PHI,
            lat_sigma: 0.06,
            
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
    noise: NoiseColorState,
    ring: RingBuffer,
    grains: GrainPool,
    envelope: Envelope,
    
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
        
        Self {
            params,
            sr,
            rng,
            w_phi,
            w_s2,
            w_pl,
            noise: NoiseColorState::default(),
            ring: RingBuffer::new(),
            grains: GrainPool::new(),
            envelope: Envelope::new(),
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
        self.noise.set_color(match self.params.noise_color {
            0 => NoiseColor::White,
            2 => NoiseColor::Brown,
            _ => NoiseColor::Pink,
        });
        self.noise.set_amount(self.params.color_amt);
        self.lattice.eps = self.params.lat_eps;
        self.lattice.gamma = self.params.lat_gamma;
        self.lattice.sigma = self.params.lat_sigma;
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
        self.noise.reset();
        self.ring.clear();
        self.grains.reset_all();
        self.lattice.reset();
        self.hawkes.reset();
        self.lat_phase = 0.0;
        self.lat_last_v = 0.0;
        self.samples_to_next = (self.sr * 0.05) as i32;
        self.gap_elapsed = 0;
        self.sample_counter = 0;
        self.lfo_wow_phase = 0.0;
        self.lfo_flut_phase = 0.0;
        self.prev_pan = 0.0;
        self.prev_itd = 0.0;
        self.prev_ild = 0.0;
        self.last_gap_samples = 0.0;
        self.last_dur_samples = 0.0;
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
    }
}

impl Engine {
    /// Internal block processing
    pub fn process_block(&mut self, out_l: &mut [f64], out_r: &mut [f64]) {
        let num_samples = out_l.len().min(out_r.len());
        
        // Pre-calculate constants
        let wow_hz = map_phi_range(0.1, 1.5, clamp01(self.params.vhs_wow));
        let flt_hz = map_phi_range(7.0, 12.0, clamp01(self.params.vhs_flutter));
        let inc_wow = wow_hz / self.sr;
        let inc_flt = flt_hz / self.sr;
        let lat_inc = clamp(self.params.lat_rate, 1.0, 2000.0) / self.sr;
        let itd_scale = self.params.itd_us * 1.0e-6 * self.sr;
        
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
                        let t = clamp01(self.params.temperature);
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
            
            // Generate and write noise to ring buffer
            let mut nz = self.noise.process(&mut self.rng);
            nz = soft_tanh(nz * 1.2);
            self.ring.write(nz);
            let wi = self.ring.get_write_index();
            
            // Schedule new grain
            self.samples_to_next -= 1;
            if self.samples_to_next <= 0 {
                self.spawn_grain(wi, vhs_mod, itd_scale);
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
                
                // Read from ring with ITD
                let itd = grain.itd + vhs_mod * 0.25 * itd_scale;
                let (mut s_l, mut s_r) = self.ring.read_stereo_itd(wi, itd);
                
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
            
            // Soft clip output
            out_l[n] = soft_tanh(y_l * OUT_DRIVE) / OUT_DRIVE;
            out_r[n] = soft_tanh(y_r * OUT_DRIVE) / OUT_DRIVE;
            
            self.sample_counter += 1;
        }
    }
    
    fn spawn_grain(&mut self, wi: usize, vhs_mod: f64, itd_scale: f64) {
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
        
        // Binaural
        let (pan_l, pan_r) = pan_equal_power(pan, self.params.width);
        grain.pan_l = pan_l;
        grain.pan_r = pan_r;
        
        let mut itd = map_itd_samples(self.sr, self.params.itd_us, pan, u2);
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
        let mut rate = clamp(self.params.rate, 0.0, MAX_EVENT_RATE_HZ);
        
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
