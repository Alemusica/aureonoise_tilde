"""
aureonoise - Offline audio analysis module

Metrics: spectral slope, stereo correlation, ITD, ILD, LUFS,
handshake rate, coherence, PLV, spectral centroid, bilateral symmetry.
"""

import numpy as np
from typing import Dict, Optional
from dataclasses import dataclass, asdict

try:
    from scipy.signal import welch, correlate
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


@dataclass
class AnalysisReport:
    """Complete analysis of a stereo audio buffer."""
    duration_sec: float
    sample_rate: float

    # Spectral
    spectral_slope: float          # dB/octave (negative = colored)
    spectral_slope_r2: float       # fit quality
    spectral_centroid_hz: float    # brightness
    spectral_spread_hz: float      # bandwidth

    # Stereo
    stereo_correlation: float      # Pearson L/R (-1 to 1)
    itd_us: float                  # measured ITD via cross-correlation
    ild_db_1k: float               # ILD at 1-2 kHz
    ild_db_4k: float               # ILD at 4-8 kHz
    bilateral_symmetry: float      # 0 = asymmetric, 1 = perfect symmetry

    # Level
    lufs_momentary: float          # ITU-R BS.1770 approximation
    peak_db: float
    rms_db: float

    # Engine metrics (optional, from Engine API)
    handshake_rate: Optional[float] = None
    coherence_mean: Optional[float] = None
    plv: Optional[float] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    def to_json(self) -> str:
        import json
        return json.dumps(self.to_dict(), indent=2)


def analyze(
    left: np.ndarray,
    right: np.ndarray,
    sample_rate: float = 44100.0,
    engine_stats: Optional[Dict] = None,
) -> AnalysisReport:
    """
    Analyze stereo audio buffer.

    Args:
        left: Left channel samples
        right: Right channel samples
        sample_rate: Sample rate in Hz
        engine_stats: Optional dict with 'handshake_rate', 'coherence', 'plv' from Engine

    Returns:
        AnalysisReport with all metrics
    """
    if not _HAS_SCIPY:
        raise ImportError("scipy required for analysis: pip install scipy")

    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    n = min(len(left), len(right))
    left, right = left[:n], right[:n]
    duration = n / sample_rate

    # --- Spectral slope (Welch PSD + linear fit in log-log) ---
    nperseg = min(2048, n // 4) if n > 256 else n
    freqs, psd_l = welch(left, fs=sample_rate, nperseg=nperseg)
    _, psd_r = welch(right, fs=sample_rate, nperseg=nperseg)
    psd_mono = (psd_l + psd_r) / 2.0

    # Fit slope in log-log space (skip DC bin)
    mask = freqs > 20.0
    if mask.sum() > 3:
        log_f = np.log10(freqs[mask])
        log_p = np.log10(psd_mono[mask] + 1e-30)
        slope, intercept = np.polyfit(log_f, log_p, 1)
        residuals = log_p - (slope * log_f + intercept)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((log_p - log_p.mean())**2)
        r2 = 1.0 - ss_res / (ss_tot + 1e-30)
        spectral_slope = slope  # dB/decade → ~slope/3.32 per octave
    else:
        spectral_slope, r2 = 0.0, 0.0

    # --- Spectral centroid and spread ---
    psd_sum = psd_mono[mask].sum()
    if psd_sum > 1e-15:
        centroid = np.sum(freqs[mask] * psd_mono[mask]) / psd_sum
        spread = np.sqrt(np.sum((freqs[mask] - centroid)**2 * psd_mono[mask]) / psd_sum)
    else:
        centroid, spread = 0.0, 0.0

    # --- Stereo correlation (Pearson) ---
    if n > 1:
        l_mean, r_mean = left.mean(), right.mean()
        l_std, r_std = left.std(), right.std()
        if l_std > 1e-12 and r_std > 1e-12:
            stereo_corr = np.mean((left - l_mean) * (right - r_mean)) / (l_std * r_std)
        else:
            stereo_corr = 1.0
    else:
        stereo_corr = 1.0

    # --- ITD measurement (cross-correlation peak, ±1.5ms window) ---
    max_lag = int(sample_rate * 0.0015)  # 1.5 ms
    if n > max_lag * 2:
        xcorr = correlate(left[:min(n, 8192)], right[:min(n, 8192)], mode='full')
        mid = len(xcorr) // 2
        lo, hi = mid - max_lag, mid + max_lag + 1
        region = xcorr[lo:hi]
        peak_idx = np.argmax(np.abs(region))
        lag_samples = peak_idx - max_lag
        itd_us = lag_samples / sample_rate * 1e6
    else:
        itd_us = 0.0

    # --- ILD per band ---
    def band_power(signal, f_lo, f_hi):
        mask_b = (freqs >= f_lo) & (freqs < f_hi)
        if mask_b.sum() == 0:
            return 1e-30
        _, psd_s = welch(signal, fs=sample_rate, nperseg=nperseg)
        return psd_s[mask_b].mean() + 1e-30

    pw_l_1k = band_power(left, 1000, 2000)
    pw_r_1k = band_power(right, 1000, 2000)
    ild_1k = 10.0 * np.log10(pw_l_1k / pw_r_1k)

    pw_l_4k = band_power(left, 4000, 8000)
    pw_r_4k = band_power(right, 4000, 8000)
    ild_4k = 10.0 * np.log10(pw_l_4k / pw_r_4k)

    # --- Bilateral symmetry ---
    # Measure how similar L and R envelopes are over time windows
    window = int(sample_rate * 0.05)  # 50ms windows
    if n > window * 2:
        n_windows = n // window
        env_l = np.array([np.sqrt(np.mean(left[i*window:(i+1)*window]**2))
                         for i in range(n_windows)])
        env_r = np.array([np.sqrt(np.mean(right[i*window:(i+1)*window]**2))
                         for i in range(n_windows)])
        env_sum = env_l + env_r + 1e-15
        symmetry = 1.0 - np.mean(np.abs(env_l - env_r) / env_sum)
    else:
        symmetry = 1.0

    # --- Level ---
    mono = (left + right) * 0.5
    rms = np.sqrt(np.mean(mono**2))
    peak = np.max(np.abs(mono))
    rms_db = 20.0 * np.log10(rms + 1e-30)
    peak_db = 20.0 * np.log10(peak + 1e-30)

    # LUFS approximation (simplified K-weighting + gating)
    # Full ITU-R BS.1770 needs K-weighting filter; this is a simplified version
    lufs = rms_db - 0.691  # approximate offset for K-weighting

    # Engine stats
    hs_rate = engine_stats.get('handshake_rate') if engine_stats else None
    coh = engine_stats.get('coherence') if engine_stats else None
    plv = engine_stats.get('plv') if engine_stats else None

    return AnalysisReport(
        duration_sec=duration,
        sample_rate=sample_rate,
        spectral_slope=spectral_slope,
        spectral_slope_r2=r2,
        spectral_centroid_hz=centroid,
        spectral_spread_hz=spread,
        stereo_correlation=stereo_corr,
        itd_us=itd_us,
        ild_db_1k=ild_1k,
        ild_db_4k=ild_4k,
        bilateral_symmetry=symmetry,
        lufs_momentary=lufs,
        peak_db=peak_db,
        rms_db=rms_db,
        handshake_rate=hs_rate,
        coherence_mean=coh,
        plv=plv,
    )
