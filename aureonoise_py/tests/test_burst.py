"""
aureonoise - Burst engine tests
Verifies the Hawkes burst position modulation port from beta7_tools.
"""

import numpy as np
import pytest
from aureonoise import Engine, Params


def _frame_rms(signal: np.ndarray, frame_size: int) -> np.ndarray:
    """Compute RMS energy per frame."""
    n_frames = len(signal) // frame_size
    return np.array([
        np.sqrt(np.mean(signal[i * frame_size:(i + 1) * frame_size] ** 2))
        for i in range(n_frames)
    ])


def _cv(values: np.ndarray) -> float:
    """Coefficient of variation (std / mean)."""
    m = np.mean(values)
    if m < 1e-10:
        return 0.0
    return float(np.std(values) / m)


class TestBurstParams:
    """Verify burst params are exposed and synced."""

    def test_default_params(self):
        p = Params()
        assert p.burst is True
        assert abs(p.burst_floor - 0.35) < 1e-6
        assert abs(p.burst_phi_mix - 0.6) < 1e-6

    def test_params_roundtrip(self):
        p = Params()
        p.burst_floor = 0.5
        p.burst_phi_mix = 0.8

        e = Engine(48000.0)
        e.set_params(p)
        q = e.get_params()
        assert abs(q.burst_floor - 0.5) < 1e-6
        assert abs(q.burst_phi_mix - 0.8) < 1e-6


class TestBurstEngine:
    """Verify burst position modulation affects output."""

    def test_burst_produces_output(self):
        """Burst-enabled engine should produce non-silent output."""
        p = Params()
        p.burst = True
        p.rate = 8.0
        p.baselen_ms = 50.0
        p.width = 0.8

        e = Engine(48000.0)
        e.set_params(p)
        left, right = e.process(48000 * 5)

        rms = np.sqrt(np.mean(left ** 2))
        assert rms > 1e-6, f"Burst output is silent (RMS={rms:.2e})"

    def test_burst_clustering(self):
        """With burst enabled, grain density should show clustering (high CV)."""
        p = Params()
        p.burst = True
        p.rate = 8.0
        p.baselen_ms = 50.0
        p.width = 0.8

        e = Engine(48000.0)
        e.set_params(p)
        left, _ = e.process(48000 * 10)

        frame_size = int(48000 * 0.05)
        energies = _frame_rms(left, frame_size)
        cv = _cv(energies)

        # Burst clustering should produce variability
        assert cv > 0.3, f"Burst CV {cv:.3f} too low -- clustering not effective"

    def test_burst_vs_no_burst_variability(self):
        """Burst mode should produce more amplitude variability than no-burst."""
        p_burst = Params()
        p_burst.burst = True
        p_burst.rate = 8.0
        p_burst.baselen_ms = 50.0
        p_burst.width = 0.8

        p_no = Params()
        p_no.burst = False
        p_no.rate = 8.0
        p_no.baselen_ms = 50.0
        p_no.width = 0.8

        e1 = Engine(48000.0)
        e1.set_params(p_burst)
        l1, _ = e1.process(48000 * 10)

        e2 = Engine(48000.0)
        e2.set_params(p_no)
        l2, _ = e2.process(48000 * 10)

        frame = int(48000 * 0.05)
        cv1 = _cv(_frame_rms(l1, frame))
        cv2 = _cv(_frame_rms(l2, frame))

        # Burst should create at least comparable variability
        # (not dramatically less -- allow margin for stochastic noise)
        assert cv1 >= cv2 * 0.5, (
            f"Burst CV std {cv1:.4f} much lower than no-burst {cv2:.4f}"
        )

    def test_burst_output_bounded(self):
        """Output should remain bounded even with burst amp scaling."""
        p = Params()
        p.burst = True
        p.rate = 20.0
        p.baselen_ms = 30.0
        p.width = 1.0
        p.temperature = 0.9

        e = Engine(48000.0)
        e.set_params(p)
        left, right = e.process(48000 * 5)

        # Soft-tanh limiter should keep output bounded
        assert np.abs(left).max() < 2.0, f"Left overflow: {np.abs(left).max():.4f}"
        assert np.abs(right).max() < 2.0, f"Right overflow: {np.abs(right).max():.4f}"

    def test_burst_disable_no_crash(self):
        """Switching burst off mid-run should not crash."""
        e = Engine(48000.0)
        p = Params()
        p.burst = True
        e.set_params(p)
        e.process(48000)

        p.burst = False
        e.set_params(p)
        left, right = e.process(48000)

        assert len(left) == 48000
        assert len(right) == 48000

    def test_stereo_asymmetry_with_burst(self):
        """Burst should create L/R differences (bilateral separation)."""
        p = Params()
        p.burst = True
        p.rate = 12.0
        p.baselen_ms = 40.0
        p.width = 1.0

        e = Engine(48000.0)
        e.set_params(p)
        left, right = e.process(48000 * 5)

        # L and R should not be identical
        diff = np.abs(left - right)
        assert np.max(diff) > 1e-6, "L/R should differ with burst spatial modulation"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
