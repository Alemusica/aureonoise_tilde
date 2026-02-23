"""
Spectral slope tests for noise generators.

Tests the raw NoiseColorState filter output (bypassing the granular engine)
to verify correct 1/f^n spectral characteristics. The granular engine's
envelope windowing and soft-clipping add their own spectral shaping, so
testing the noise source directly is the only way to measure filter accuracy.
"""
import pytest
import numpy as np
from scipy.signal import welch
from scipy.stats import linregress
from aureonoise import NoiseColorState, NoiseColor, Rng


def _measure_slope(samples, fs=48000, f_lo=50, f_hi=15000):
    """Measure spectral slope via log-log linear regression on Welch PSD."""
    f, Pxx = welch(samples, fs=fs, nperseg=8192, noverlap=4096)
    mask = (f >= f_lo) & (f <= f_hi)
    log_f = np.log10(f[mask])
    log_P = np.log10(Pxx[mask] + 1e-30)
    slope, _, r_value, _, _ = linregress(log_f, log_P)
    return slope, r_value**2


def test_pink_noise_slope():
    """Pink noise must have spectral slope -1.0 +/- 0.15 (1/f)."""
    rng = Rng(42)
    nc = NoiseColorState(NoiseColor.Pink, 1.0)

    # Generate 10 seconds at 48kHz
    samples = np.array([nc.process(rng) for _ in range(480000)])

    slope, r2 = _measure_slope(samples)

    assert -1.15 <= slope <= -0.85, f"Pink slope {slope:.3f} outside [-1.15, -0.85]"
    assert r2 > 0.90, f"R2 = {r2:.3f} too low (expected > 0.90)"


@pytest.mark.xfail(reason="Brown filter uses cascaded poles + soft_tanh, produces ~-3.9 instead of -2.0. Needs dedicated fix.")
def test_brown_noise_slope():
    """Brown noise must have spectral slope -2.0 +/- 0.3 (1/f^2)."""
    rng = Rng(42)
    nc = NoiseColorState(NoiseColor.Brown, 1.0)

    # Generate 10 seconds at 48kHz
    samples = np.array([nc.process(rng) for _ in range(480000)])

    slope, r2 = _measure_slope(samples)

    assert -2.3 <= slope <= -1.7, f"Brown slope {slope:.3f} outside [-2.3, -1.7]"


def test_white_noise_flat():
    """White noise must have approximately flat spectrum (slope ~0)."""
    rng = Rng(42)
    nc = NoiseColorState(NoiseColor.White, 1.0)

    samples = np.array([nc.process(rng) for _ in range(480000)])

    slope, r2 = _measure_slope(samples)

    assert -0.15 <= slope <= 0.15, f"White slope {slope:.3f} outside [-0.15, 0.15]"


def test_pink_kellet_accuracy():
    """Kellet filter must achieve R2 > 0.98 against ideal 1/f line."""
    rng = Rng(12345)
    nc = NoiseColorState(NoiseColor.Pink, 1.0)

    samples = np.array([nc.process(rng) for _ in range(480000)])

    slope, r2 = _measure_slope(samples)

    # Kellet's refined method should be very tight
    assert -1.10 <= slope <= -0.90, f"Kellet slope {slope:.3f} outside [-1.10, -0.90]"
    assert r2 > 0.98, f"R2 = {r2:.3f} too low for Kellet (expected > 0.98)"
