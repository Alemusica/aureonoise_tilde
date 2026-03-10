"""
aureonoise - Audio analysis module

Measures signal quality and therapeutic effectiveness.
Metrics: spectral slope, stereo correlation, ITD, ILD, LUFS,
handshake rate, coherence, PLV, spectral centroid, bilateral symmetry.

All spectral slopes in dB/octave (not dB/decade).
Pink = -1 dB/oct, Brown = -2 dB/oct, White = 0 dB/oct.
"""

import numpy as np
from typing import Dict, Optional, List
from dataclasses import dataclass, asdict, field

try:
    from scipy.signal import welch, correlate
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


# Therapeutic quality thresholds
QUALITY_THRESHOLDS = {
    "stereo_correlation": {"good": (-0.3, 0.6), "warn": "Outside range: mono or phase issues"},
    "itd_us": {"good": (0, 800), "warn": "ITD > 800 us: unnatural spatial cue"},
    "bilateral_symmetry": {"good": (0.3, 0.95), "warn": "Outside range: asymmetric L/R"},
    "peak_db": {"good": (-40, -3), "warn": "Peak too hot or too quiet"},
    "spectral_slope_r2": {"good": (0.7, 1.0), "warn": "Poor spectral fit: noisy or aliased"},
}


@dataclass
class AnalysisReport:
    """Complete analysis of a stereo audio buffer."""
    duration_sec: float
    sample_rate: float

    # Spectral
    spectral_slope: float          # dB/octave (0=white, -1=pink, -2=brown)
    spectral_slope_r2: float       # fit quality [0,1]
    spectral_centroid_hz: float    # brightness
    spectral_spread_hz: float      # bandwidth

    # Stereo
    stereo_correlation: float      # Pearson L/R (-1 to 1)
    itd_us: float                  # measured ITD via cross-correlation (microseconds)
    ild_db_1k: float               # ILD at 1-2 kHz (dB)
    ild_db_4k: float               # ILD at 4-8 kHz (dB)
    bilateral_symmetry: float      # 0 = asymmetric, 1 = perfect symmetry

    # Level
    lufs_momentary: float          # simplified ITU-R BS.1770
    peak_db: float
    rms_db: float

    # Engine metrics (optional, from Engine API)
    handshake_rate: Optional[float] = None
    coherence_mean: Optional[float] = None
    plv: Optional[float] = None

    # Quality warnings
    warnings: List[str] = field(default_factory=list)
    verdict: str = ""              # "OK", "WARN", "FAIL"

    def to_dict(self) -> Dict:
        return asdict(self)

    def to_json(self) -> str:
        import json
        return json.dumps(self.to_dict(), indent=2)

    def summary_lines(self) -> List[str]:
        """Human-readable summary for GUI display."""
        lines = [
            f"Slope: {self.spectral_slope:.2f} dB/oct (R2={self.spectral_slope_r2:.2f})",
            f"Centroid: {self.spectral_centroid_hz:.0f} Hz  Spread: {self.spectral_spread_hz:.0f} Hz",
            f"Stereo corr: {self.stereo_correlation:.3f}",
            f"ITD: {self.itd_us:.1f} us",
            f"ILD: 1k={self.ild_db_1k:.1f} dB  4k={self.ild_db_4k:.1f} dB",
            f"Bilateral sym: {self.bilateral_symmetry:.3f}",
            f"Level: {self.rms_db:.1f} dBFS (peak {self.peak_db:.1f})",
        ]
        if self.handshake_rate is not None:
            lines.append(f"Handshakes: {self.handshake_rate:.3f}")
        if self.coherence_mean is not None:
            lines.append(f"Coherence: {self.coherence_mean:.3f}")
        if self.plv is not None:
            lines.append(f"PLV: {self.plv:.3f}")
        return lines


def analyze(
    left: np.ndarray,
    right: np.ndarray,
    sample_rate: float = 44100.0,
    engine_stats: Optional[Dict] = None,
) -> AnalysisReport:
    """
    Analyze stereo audio buffer.

    Args:
        left: Left channel samples (float64 or float32)
        right: Right channel samples
        sample_rate: Sample rate in Hz
        engine_stats: Optional dict with 'handshake_rate', 'coherence', 'plv' from Engine

    Returns:
        AnalysisReport with all metrics and quality verdict
    """
    if not _HAS_SCIPY:
        raise ImportError("scipy required for analysis: pip install scipy")

    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    n = min(len(left), len(right))
    if n < 64:
        return _empty_report(n / sample_rate, sample_rate)

    left, right = left[:n], right[:n]
    duration = n / sample_rate

    # --- Welch PSD (computed once, reused for all spectral metrics) ---
    nperseg = min(4096, n // 2) if n > 512 else n
    freqs, psd_l = welch(left, fs=sample_rate, nperseg=nperseg)
    _, psd_r = welch(right, fs=sample_rate, nperseg=nperseg)
    psd_mono = (psd_l + psd_r) / 2.0

    # --- Spectral slope in dB/OCTAVE ---
    # Fit in log2(f) space so slope reads directly as dB/octave.
    # Pink = -1 dB/oct, Brown = -2, White = 0
    mask = freqs > 30.0
    if mask.sum() > 5:
        log2_f = np.log2(freqs[mask])
        db_p = 10.0 * np.log10(psd_mono[mask] + 1e-30)
        slope, intercept = np.polyfit(log2_f, db_p, 1)
        residuals = db_p - (slope * log2_f + intercept)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((db_p - db_p.mean()) ** 2)
        r2 = float(1.0 - ss_res / (ss_tot + 1e-30))
        spectral_slope = float(slope)
    else:
        spectral_slope, r2 = 0.0, 0.0

    # --- Spectral centroid and spread ---
    if mask.sum() > 0:
        f_mask = freqs[mask]
        p_mask = psd_mono[mask]
        psd_sum = p_mask.sum()
        if psd_sum > 1e-15:
            centroid = float(np.sum(f_mask * p_mask) / psd_sum)
            spread = float(np.sqrt(np.sum((f_mask - centroid) ** 2 * p_mask) / psd_sum))
        else:
            centroid, spread = 0.0, 0.0
    else:
        centroid, spread = 0.0, 0.0

    # --- Stereo correlation (Pearson) ---
    if n > 1:
        l_std, r_std = left.std(), right.std()
        if l_std > 1e-12 and r_std > 1e-12:
            stereo_corr = float(np.corrcoef(left, right)[0, 1])
        else:
            stereo_corr = 1.0
    else:
        stereo_corr = 1.0

    # --- ITD measurement (cross-correlation peak, +/-1.5ms window) ---
    max_lag = int(sample_rate * 0.0015)
    if n > max_lag * 2 + 64:
        chunk = min(n, 16384)
        xcorr = correlate(left[:chunk], right[:chunk], mode='full')
        mid = len(xcorr) // 2
        lo, hi = mid - max_lag, mid + max_lag + 1
        region = xcorr[lo:hi]
        peak_idx = int(np.argmax(np.abs(region)))
        lag_samples = peak_idx - max_lag
        itd_us = float(lag_samples / sample_rate * 1e6)
    else:
        itd_us = 0.0

    # --- ILD per band (reuse pre-computed PSD, no extra Welch calls) ---
    ild_1k = _band_ild(freqs, psd_l, psd_r, 1000, 2000)
    ild_4k = _band_ild(freqs, psd_l, psd_r, 4000, 8000)

    # --- Bilateral symmetry ---
    window = int(sample_rate * 0.05)  # 50ms windows
    if n > window * 4:
        n_windows = n // window
        env_l = np.array([np.sqrt(np.mean(left[i * window:(i + 1) * window] ** 2))
                          for i in range(n_windows)])
        env_r = np.array([np.sqrt(np.mean(right[i * window:(i + 1) * window] ** 2))
                          for i in range(n_windows)])
        env_sum = env_l + env_r + 1e-15
        symmetry = float(1.0 - np.mean(np.abs(env_l - env_r) / env_sum))
    else:
        symmetry = 1.0

    # --- Level ---
    mono = (left + right) * 0.5
    rms = float(np.sqrt(np.mean(mono ** 2)))
    peak = float(np.max(np.abs(mono)))
    rms_db = 20.0 * np.log10(rms + 1e-30)
    peak_db = 20.0 * np.log10(peak + 1e-30)
    lufs = rms_db - 0.691

    # Engine stats
    hs_rate = engine_stats.get('handshake_rate') if engine_stats else None
    coh = engine_stats.get('coherence') if engine_stats else None
    plv_val = engine_stats.get('plv') if engine_stats else None

    # --- Quality verdict ---
    warnings = _check_quality(stereo_corr, abs(itd_us), symmetry, peak_db, r2)
    verdict = "OK" if not warnings else ("FAIL" if len(warnings) > 2 else "WARN")

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
        plv=plv_val,
        warnings=warnings,
        verdict=verdict,
    )


def _band_ild(freqs, psd_l, psd_r, f_lo, f_hi):
    """ILD for a frequency band using pre-computed PSD arrays."""
    band = (freqs >= f_lo) & (freqs < f_hi)
    if band.sum() == 0:
        return 0.0
    pw_l = psd_l[band].mean() + 1e-30
    pw_r = psd_r[band].mean() + 1e-30
    return float(10.0 * np.log10(pw_l / pw_r))


def _check_quality(stereo_corr, itd_abs, symmetry, peak_db, slope_r2):
    """Check signal against therapeutic quality thresholds."""
    warns = []
    if stereo_corr > 0.95:
        warns.append("Stereo too correlated (mono-like)")
    elif stereo_corr < -0.5:
        warns.append("Stereo anti-correlated: phase cancellation")
    if itd_abs > 900:
        warns.append(f"ITD {itd_abs:.0f} us exceeds physiological range")
    if symmetry < 0.2:
        warns.append(f"Bilateral symmetry {symmetry:.2f}: heavily asymmetric")
    if peak_db > -1:
        warns.append(f"Peak {peak_db:.1f} dBFS: clipping risk")
    elif peak_db < -50:
        warns.append(f"Peak {peak_db:.1f} dBFS: signal too quiet")
    if slope_r2 < 0.5:
        warns.append(f"Spectral fit R2={slope_r2:.2f}: noisy spectrum")
    return warns


def _empty_report(duration, sr):
    """Return zeroed report for buffers too short to analyze."""
    return AnalysisReport(
        duration_sec=duration, sample_rate=sr,
        spectral_slope=0, spectral_slope_r2=0,
        spectral_centroid_hz=0, spectral_spread_hz=0,
        stereo_correlation=0, itd_us=0,
        ild_db_1k=0, ild_db_4k=0,
        bilateral_symmetry=0,
        lufs_momentary=-80, peak_db=-80, rms_db=-80,
        warnings=["Buffer too short for analysis"],
        verdict="FAIL",
    )
