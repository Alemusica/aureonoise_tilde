"""
aureonoise - Test suite
"""

import pytest
import numpy as np
from aureonoise import (
    Engine, Params, Rng, Weyl,
    PHI, INV_PHI, INV_PHI_SQ, INV_PHI_CU,
    MAX_GRAINS, RING_SIZE,
)
from aureonoise.presets import PresetBank, Preset


class TestConstants:
    """Test φ-based constants."""
    
    def test_phi_values(self):
        assert abs(PHI - 1.6180339887498948) < 1e-10
        assert abs(INV_PHI - 0.6180339887498948) < 1e-10
        assert abs(INV_PHI_SQ - 0.3819660112501051) < 1e-10
        assert abs(INV_PHI_CU - 0.2360679774997896) < 1e-10
    
    def test_phi_identity(self):
        # φ - 1 = 1/φ
        assert abs((PHI - 1) - INV_PHI) < 1e-10
        # 1/φ - 1/φ² = 1/φ³
        assert abs((INV_PHI - INV_PHI_SQ) - INV_PHI_CU) < 1e-10
    
    def test_system_constants(self):
        assert MAX_GRAINS == 32
        assert RING_SIZE == 131072


class TestRng:
    """Test random number generator."""
    
    def test_creation(self):
        rng = Rng(12345)
        assert rng is not None
    
    def test_uni01_range(self):
        rng = Rng(42)
        for _ in range(1000):
            v = rng.uni01()
            assert 0.0 <= v <= 1.0
    
    def test_reproducibility(self):
        rng1 = Rng(42)
        rng2 = Rng(42)
        for _ in range(100):
            assert rng1.uni01() == rng2.uni01()


class TestWeyl:
    """Test Weyl sequence generators."""
    
    def test_phi_sequence(self):
        w = Weyl.phi(0.0)
        values = [w.next() for _ in range(1000)]
        assert all(0.0 <= v <= 1.0 for v in values)
    
    def test_coverage(self):
        w = Weyl.phi(0.0)
        bins = [0] * 10
        for _ in range(10000):
            v = w.next()
            bins[int(v * 9.999)] += 1
        # All bins should have roughly equal counts
        avg = sum(bins) / len(bins)
        for b in bins:
            assert abs(b - avg) < avg * 0.2


class TestEngine:
    """Test DSP engine."""
    
    def test_creation(self):
        engine = Engine(44100.0)
        assert engine.sample_rate() == 44100.0
    
    def test_params(self):
        engine = Engine(44100.0)
        params = engine.get_params()
        assert params.rate == 8.0
        assert params.baselen_ms == 120.0
    
    def test_set_params(self):
        engine = Engine(44100.0)
        params = Params()
        params.rate = 15.0
        params.glitch_mix = 0.8
        engine.set_params(params)
        
        new_params = engine.get_params()
        assert new_params.rate == 15.0
        assert new_params.glitch_mix == 0.8
    
    def test_process(self):
        engine = Engine(44100.0)
        left, right = engine.process(1024)
        
        assert len(left) == 1024
        assert len(right) == 1024
        assert isinstance(left, np.ndarray)
        assert isinstance(right, np.ndarray)
    
    def test_output_range(self):
        engine = Engine(44100.0)
        params = Params()
        params.rate = 20.0
        engine.set_params(params)
        
        left, right = engine.process(44100)
        
        # Output should be bounded
        assert np.abs(left).max() < 2.0
        assert np.abs(right).max() < 2.0
    
    def test_reset(self):
        engine = Engine(44100.0)
        engine.process(10000)
        engine.reset()
        
        # Should be able to process again
        left, right = engine.process(1000)
        assert len(left) == 1000
    
    def test_different_sample_rates(self):
        for sr in [22050.0, 44100.0, 48000.0, 96000.0]:
            engine = Engine(sr)
            assert engine.sample_rate() == sr
            left, right = engine.process(1000)
            assert len(left) == 1000


class TestPresets:
    """Test preset system."""
    
    def test_preset_bank(self):
        bank = PresetBank()
        names = bank.list()
        
        assert "Default" in names
        assert "Digital Storm" in names
        assert len(names) >= 5
    
    def test_get_preset(self):
        bank = PresetBank()
        preset = bank.get("Default")
        
        assert preset is not None
        assert preset.name == "Default"
        assert len(preset.params) > 0
    
    def test_to_params(self):
        bank = PresetBank()
        preset = bank.get("Digital Storm")
        params = preset.to_params()
        
        assert params.rate == 20.0
        assert params.glitch_mix == 0.9
    
    def test_from_params(self):
        params = Params()
        params.rate = 25.0
        params.temperature = 0.9
        
        preset = Preset.from_params("Custom", params, "Test preset")
        
        assert preset.name == "Custom"
        assert preset.params["rate"] == 25.0
        assert preset.params["temperature"] == 0.9
    
    def test_json_roundtrip(self):
        bank = PresetBank()
        original = bank.get("Calm Rain")
        
        json_str = original.to_json()
        loaded = Preset.from_json(json_str)
        
        assert loaded.name == original.name
        assert loaded.params["rate"] == original.params["rate"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
