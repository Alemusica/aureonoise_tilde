"""Smoke test: app instantiation and DPG context creation."""
import pytest


def test_app_instantiation():
    """AureonoiseApp can be created without crashing."""
    from aureonoise.app import AureonoiseApp
    app = AureonoiseApp(sample_rate=44100.0, block_size=512)
    assert app.audio is not None
    assert app.running is False


def test_imports():
    """All app.py imports resolve."""
    from aureonoise.app import COLORS, NOISE_MODES, MODAL_PRESETS, main
    assert len(COLORS) > 0
    assert len(NOISE_MODES) == 6
    assert len(MODAL_PRESETS) == 4
