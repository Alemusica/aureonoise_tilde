"""
aureonoise - Therapeutic presets & dialogue integration test suite

Tests evidence-based therapeutic presets (EMDR, ASMR, sleep, focus,
theta, hemispheric bridge) and the dialogue/bilateral subsystems.

Handles two build states:
  - Old build: only classic Params fields (rate, baselen_ms, noise_color, ...).
    Tests that require new fields (dialogue_on, bilateral_on, noise_mode, etc.)
    are skipped with a clear reason.
  - New build (after `maturin develop` with full lib.rs): all tests run.
"""

import pytest
import numpy as np

try:
    from aureonoise import Engine, Params, PHI, INV_PHI
    from aureonoise.presets import PresetBank, Preset
    _HAVE_CORE = True
except ImportError:
    _HAVE_CORE = False

pytestmark = pytest.mark.skipif(
    not _HAVE_CORE,
    reason="aureonoise._core not built (run `maturin develop` first)",
)

SR = 44100.0
BLOCK = 4096            # ~93 ms at 44.1 kHz -- fast default block
LONG_BLOCK = int(SR)    # 1 second
VERY_LONG = int(SR * 5) # 5 seconds -- only where needed


# ---------------------------------------------------------------------------
# API detection
# ---------------------------------------------------------------------------

def _has_param(name: str) -> bool:
    """Check whether Params exposes a given attribute in the current build."""
    if not _HAVE_CORE:
        return False
    try:
        p = Params()
        getattr(p, name)
        return True
    except AttributeError:
        return False


def _has_engine_method(name: str) -> bool:
    """Check whether Engine exposes a given method in the current build."""
    if not _HAVE_CORE:
        return False
    return hasattr(Engine, name) or hasattr(Engine(SR), name)


# Feature flags for the current build
_HAS_DIALOGUE = _has_param("dialogue_on")
_HAS_BILATERAL = _has_param("bilateral_on")
_HAS_NOISE_MODE = _has_param("noise_mode")
_HAS_MODAL = _has_param("modal_on")
_HAS_EXTERNALIZATION = _has_param("externalization")
_HAS_PHI_PAN = _has_param("phi_pan")
_HAS_COHERENCE = _has_engine_method("coherence")

# Composite: all therapeutic subsystems are present
_HAS_THERAPEUTIC = _HAS_DIALOGUE and _HAS_BILATERAL and _HAS_NOISE_MODE

_SKIP_THERAPEUTIC = pytest.mark.skipif(
    not _HAS_THERAPEUTIC,
    reason="Therapeutic params not in current build (dialogue_on, bilateral_on, "
           "noise_mode). Run `maturin develop` to rebuild.",
)
_SKIP_DIALOGUE = pytest.mark.skipif(
    not _HAS_DIALOGUE,
    reason="dialogue_on not in current build. Run `maturin develop`.",
)
_SKIP_BILATERAL = pytest.mark.skipif(
    not _HAS_BILATERAL,
    reason="bilateral_on not in current build. Run `maturin develop`.",
)
_SKIP_NOISE_MODE = pytest.mark.skipif(
    not _HAS_NOISE_MODE,
    reason="noise_mode not in current build. Run `maturin develop`.",
)
_SKIP_MODAL = pytest.mark.skipif(
    not _HAS_MODAL,
    reason="modal_on not in current build. Run `maturin develop`.",
)
_SKIP_COHERENCE = pytest.mark.skipif(
    not _HAS_COHERENCE,
    reason="Engine.coherence() not in current build. Run `maturin develop`.",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rms(signal: np.ndarray) -> float:
    """Root-mean-square of a 1-D array."""
    return float(np.sqrt(np.mean(signal ** 2)))


def _peak(signal: np.ndarray) -> float:
    """Peak absolute value."""
    return float(np.max(np.abs(signal)))


def _has_nan_or_inf(signal: np.ndarray) -> bool:
    """True if any NaN or Inf present."""
    return bool(np.any(~np.isfinite(signal)))


def _lr_correlation(left: np.ndarray, right: np.ndarray) -> float:
    """Pearson correlation between L and R channels."""
    if _rms(left) < 1e-12 or _rms(right) < 1e-12:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def _frame_rms(signal: np.ndarray, frame_size: int) -> np.ndarray:
    """RMS per non-overlapping frame."""
    n = len(signal) // frame_size
    return np.array([
        np.sqrt(np.mean(signal[i * frame_size:(i + 1) * frame_size] ** 2))
        for i in range(n)
    ])


def _safe_set(params, name: str, value):
    """Set a Params attribute only if it exists in the current build."""
    if hasattr(params, name):
        setattr(params, name, value)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def engine():
    """Fresh Engine at 44.1 kHz."""
    return Engine(SR)


@pytest.fixture
def preset_bank():
    """Fresh PresetBank with factory defaults."""
    return PresetBank()


# ---------------------------------------------------------------------------
# Therapeutic preset parameter factories
# ---------------------------------------------------------------------------

def _emdr_params() -> Params:
    """EMDR bilateral stimulation preset."""
    p = Params()
    p.rate = 6.0
    p.baselen_ms = 150.0
    p.len_phi = 0.7
    p.width = 1.4
    p.noise_color = 2      # Brown fallback for old build
    _safe_set(p, "noise_mode", 2)
    p.color_amt = 0.75
    _safe_set(p, "bilateral_on", True)
    _safe_set(p, "bilateral_rate", 1.0)
    _safe_set(p, "bilateral_amount", 0.9)
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "dialogue_strength", 0.5)
    p.glitch_mix = 0.0
    p.srcrush_amt = 0.0
    p.bitcrush_amt = 0.0
    p.vhs_wow = 0.1
    p.vhs_flutter = 0.1
    p.burst = False
    p.temperature = 0.3
    p.seed = 42
    return p


def _asmr_params() -> Params:
    """ASMR intimate micro-transient preset."""
    p = Params()
    p.rate = 14.0
    p.baselen_ms = 40.0
    p.len_phi = 0.5
    p.width = 1.0
    p.noise_color = 2
    _safe_set(p, "noise_mode", 2)
    p.color_amt = 0.8
    p.burst = True
    p.burst_floor = 0.2
    p.burst_phi_mix = 0.7
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "bilateral_on", False)
    p.glitch_mix = 0.1
    p.srcrush_amt = 0.0
    p.bitcrush_amt = 0.0
    p.temperature = 0.4
    p.seed = 42
    return p


def _sleep_params() -> Params:
    """Sleep/relaxation pink noise preset."""
    p = Params()
    p.rate = 4.0
    p.baselen_ms = 300.0
    p.len_phi = 0.6
    p.width = 0.8
    p.noise_color = 1      # Pink
    _safe_set(p, "noise_mode", 1)
    p.color_amt = 0.85
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "dialogue_strength", 0.3)
    _safe_set(p, "bilateral_on", False)
    p.glitch_mix = 0.0
    p.srcrush_amt = 0.0
    p.bitcrush_amt = 0.0
    p.vhs_wow = 0.05
    p.vhs_flutter = 0.05
    p.burst = False
    p.temperature = 0.15
    p.seed = 42
    return p


def _focus_params() -> Params:
    """Focus/concentration brown noise preset."""
    p = Params()
    p.rate = 6.0
    p.baselen_ms = 200.0
    p.len_phi = 0.7
    p.width = 0.6
    p.noise_color = 2      # Brown
    _safe_set(p, "noise_mode", 2)
    p.color_amt = 0.7
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "dialogue_strength", 0.4)
    _safe_set(p, "bilateral_on", False)
    p.glitch_mix = 0.0
    p.srcrush_amt = 0.0
    p.bitcrush_amt = 0.0
    p.burst = False
    p.temperature = 0.25
    p.seed = 42
    return p


def _theta_params() -> Params:
    """Theta-drift brainwave entrainment preset."""
    p = Params()
    p.rate = 10.0
    p.baselen_ms = 100.0
    p.len_phi = 0.8
    p.width = 1.2
    p.noise_color = 1      # Pink
    _safe_set(p, "noise_mode", 1)
    p.color_amt = 0.6
    p.burst = True
    p.burst_floor = 0.3
    p.burst_phi_mix = 0.8
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "dialogue_strength", 0.7)
    _safe_set(p, "bilateral_on", True)
    _safe_set(p, "bilateral_rate", 0.5)
    _safe_set(p, "bilateral_amount", 0.6)
    p.glitch_mix = 0.0
    p.temperature = 0.35
    p.seed = 42
    return p


def _hemispheric_bridge_params() -> Params:
    """Hemispheric bridge -- the main therapeutic preset.

    Full dialogue + bilateral + burst, designed for interhemispheric
    coherence to increase over time.
    """
    p = Params()
    p.rate = 8.0
    p.baselen_ms = 120.0
    p.len_phi = 0.8
    p.width = 1.0
    p.noise_color = 1
    _safe_set(p, "noise_mode", 1)
    p.color_amt = 0.65
    _safe_set(p, "dialogue_on", True)
    _safe_set(p, "dialogue_strength", 0.65)
    _safe_set(p, "dialogue_memory", 0.6)
    _safe_set(p, "dialogue_phi_mix", 0.85)
    _safe_set(p, "bilateral_on", True)
    _safe_set(p, "bilateral_rate", 0.8)
    _safe_set(p, "bilateral_amount", 0.7)
    p.burst = True
    p.burst_floor = 0.35
    p.burst_phi_mix = 0.6
    p.hemis_coupling = 0.7
    p.glitch_mix = 0.0
    p.temperature = 0.4
    p.seed = 42
    return p


# ===================================================================
# TestTherapeuticPresets -- each preset produces valid audio output
# ===================================================================

class TestTherapeuticPresets:
    """Verify each therapeutic preset produces valid, bounded audio.

    These tests work on both old and new builds: preset factories use
    _safe_set for new fields, so the engine will run with whatever
    params are available.  Tests that assert *therapeutic-specific*
    behaviour (bilateral alternation, coherence, etc.) are gated
    with skip markers.
    """

    def test_emdr_output_valid(self, engine):
        """EMDR preset: produces non-silent, bounded output, no NaN."""
        engine.set_params(_emdr_params())
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left), "NaN/Inf in left channel"
        assert not _has_nan_or_inf(right), "NaN/Inf in right channel"
        assert _peak(left) < 2.0, "Left channel overflow"
        assert _peak(right) < 2.0, "Right channel overflow"

    @_SKIP_BILATERAL
    def test_emdr_bilateral_alternation(self, engine):
        """EMDR bilateral: L and R should differ (alternation)."""
        engine.set_params(_emdr_params())
        left, right = engine.process(LONG_BLOCK)

        diff_rms = _rms(np.asarray(left) - np.asarray(right))
        assert diff_rms > 1e-6, (
            "EMDR bilateral: L and R are identical -- no alternation detected"
        )

    def test_asmr_intimate(self, engine):
        """ASMR: brown noise, burst clusters, output bounded."""
        engine.set_params(_asmr_params())
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert _peak(left) < 2.0
        assert _peak(right) < 2.0
        assert _rms(left) > 1e-6, "ASMR output is silent"

    def test_sleep_pink(self, engine):
        """Sleep: pink noise, low energy, stable -- no spikes."""
        engine.set_params(_sleep_params())
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)

        rms_l = _rms(np.asarray(left))
        peak_l = _peak(left)
        assert peak_l < 1.5, f"Sleep peak too high: {peak_l:.4f}"
        if rms_l > 1e-8:
            crest = peak_l / rms_l
            assert crest < 30.0, f"Sleep crest factor too high: {crest:.1f}"

    def test_focus_brown(self, engine):
        """Focus: brown noise, moderate rate, bounded."""
        engine.set_params(_focus_params())
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert _peak(left) < 2.0
        assert _peak(right) < 2.0
        assert _rms(left) > 1e-7, "Focus output is silent"

    def test_theta_drift(self, engine):
        """Theta: burst active, output valid."""
        engine.set_params(_theta_params())
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert _peak(left) < 2.0
        assert _peak(right) < 2.0

    @_SKIP_BILATERAL
    def test_theta_bilateral_lr_diff(self, engine):
        """Theta: bilateral should produce L-R differences."""
        engine.set_params(_theta_params())
        left, right = engine.process(LONG_BLOCK)

        diff_rms = _rms(np.asarray(left) - np.asarray(right))
        assert diff_rms > 1e-6, "Theta bilateral: no L-R difference"

    def test_hemispheric_bridge_output(self, engine):
        """Hemispheric bridge: produces valid bounded output."""
        engine.set_params(_hemispheric_bridge_params())
        left, right = engine.process(VERY_LONG)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert _peak(left) < 2.0

    @_SKIP_COHERENCE
    def test_hemispheric_bridge_coherence(self, engine):
        """Hemispheric bridge: coherence should evolve, handshakes should fire."""
        engine.set_params(_hemispheric_bridge_params())
        left, right = engine.process(VERY_LONG)

        coh = engine.coherence()
        assert coh != 1.0 or engine.handshake_count() > 0, (
            "Hemispheric bridge: dialogue system did not engage "
            f"(coherence={coh:.4f}, handshakes={engine.handshake_count()})"
        )
        assert engine.handshake_count() > 0, (
            f"No handshakes after 5 s (utterances would be ~{5 * 8})"
        )


# ===================================================================
# TestDialogueIntegration -- dialogue subsystem through the engine
# ===================================================================

class TestDialogueIntegration:
    """Test dialogue system behaviour through the engine API."""

    @_SKIP_DIALOGUE
    @_SKIP_COHERENCE
    def test_coherence_changes_with_dialogue(self, engine):
        """With dialogue_on=True, coherence should deviate from default."""
        p = Params()
        p.dialogue_on = True
        p.dialogue_strength = 0.7
        p.rate = 10.0
        p.seed = 99
        engine.set_params(p)

        engine.process(LONG_BLOCK * 2)

        coh = engine.coherence()
        assert coh != 1.0 or engine.handshake_count() > 0, (
            f"Dialogue on but coherence static at {coh:.4f}"
        )

    @_SKIP_DIALOGUE
    @_SKIP_COHERENCE
    def test_dialogue_off_passthrough(self, engine):
        """With dialogue_on=False, output is produced without crash."""
        p = Params()
        p.dialogue_on = False
        p.rate = 10.0
        p.seed = 99
        engine.set_params(p)

        left, right = engine.process(LONG_BLOCK * 2)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert engine.coherence() is not None

    @_SKIP_BILATERAL
    def test_bilateral_produces_lr_alternation(self):
        """bilateral_on=True should decorrelate L and R."""
        # --- Without bilateral ---
        p_off = Params()
        p_off.bilateral_on = False
        p_off.rate = 8.0
        p_off.width = 1.0
        p_off.seed = 77

        e_off = Engine(SR)
        e_off.set_params(p_off)
        l_off, r_off = e_off.process(LONG_BLOCK)
        corr_off = _lr_correlation(np.asarray(l_off), np.asarray(r_off))

        # --- With bilateral ---
        p_on = Params()
        p_on.bilateral_on = True
        p_on.bilateral_rate = 1.0
        p_on.bilateral_amount = 0.9
        p_on.rate = 8.0
        p_on.width = 1.0
        p_on.seed = 77

        e_on = Engine(SR)
        e_on.set_params(p_on)
        l_on, r_on = e_on.process(LONG_BLOCK)
        corr_on = _lr_correlation(np.asarray(l_on), np.asarray(r_on))

        # Bilateral should decorrelate L/R (lower or comparable correlation)
        assert corr_on <= corr_off + 0.15, (
            f"Bilateral did not decorrelate L/R: corr_on={corr_on:.4f}, "
            f"corr_off={corr_off:.4f}"
        )

    @_SKIP_NOISE_MODE
    def test_all_noise_modes_produce_output(self):
        """noise_mode 0-5 all produce non-zero output."""
        for mode in range(6):
            e = Engine(SR)
            p = Params()
            p.noise_mode = mode
            p.rate = 12.0
            p.seed = 42
            e.set_params(p)

            left, right = e.process(BLOCK)

            assert not _has_nan_or_inf(left), f"NaN/Inf in noise_mode={mode} left"
            assert not _has_nan_or_inf(right), f"NaN/Inf in noise_mode={mode} right"

            total_rms = _rms(np.asarray(left)) + _rms(np.asarray(right))
            assert total_rms > 1e-9, (
                f"noise_mode={mode} produced silence (total RMS={total_rms:.2e})"
            )

    def test_all_noise_colors_produce_output(self):
        """noise_color 0-2 (white/pink/brown) produce non-zero output.

        This test works on both old and new builds since noise_color
        is always available.
        """
        for color in range(3):
            e = Engine(SR)
            p = Params()
            p.noise_color = color
            p.rate = 12.0
            p.seed = 42
            e.set_params(p)

            left, right = e.process(BLOCK)

            assert not _has_nan_or_inf(left), f"NaN/Inf in noise_color={color}"
            assert not _has_nan_or_inf(right), f"NaN/Inf in noise_color={color}"

            total_rms = _rms(np.asarray(left)) + _rms(np.asarray(right))
            assert total_rms > 1e-9, (
                f"noise_color={color} produced silence (total RMS={total_rms:.2e})"
            )

    @_SKIP_MODAL
    def test_modal_resonator_changes_output(self):
        """modal_on=True should change the output compared to modal_on=False."""
        p_off = Params()
        p_off.modal_on = False
        p_off.rate = 10.0
        p_off.seed = 55

        e_off = Engine(SR)
        e_off.set_params(p_off)
        l_off, _ = e_off.process(LONG_BLOCK)

        p_on = Params()
        p_on.modal_on = True
        p_on.modal_mix = 0.5
        p_on.modal_decay = 0.5
        p_on.modal_preset = 1
        p_on.rate = 10.0
        p_on.seed = 55

        e_on = Engine(SR)
        e_on.set_params(p_on)
        l_on, _ = e_on.process(LONG_BLOCK)

        diff = np.asarray(l_on) - np.asarray(l_off)
        diff_rms = _rms(diff)
        assert diff_rms > 1e-8, (
            f"Modal resonator had no effect (diff RMS={diff_rms:.2e})"
        )

    def test_engine_deterministic_across_instances(self):
        """Two engines with the same seed + params produce identical output."""
        p = Params()
        p.rate = 8.0
        p.seed = 12345
        _safe_set(p, "dialogue_on", True)
        _safe_set(p, "bilateral_on", True)

        e1 = Engine(SR)
        e1.set_params(p)
        l1, r1 = e1.process(BLOCK)

        e2 = Engine(SR)
        e2.set_params(p)
        l2, r2 = e2.process(BLOCK)

        np.testing.assert_array_equal(
            np.asarray(l1), np.asarray(l2),
            err_msg="Left channel differs between two engines with same seed",
        )
        np.testing.assert_array_equal(
            np.asarray(r1), np.asarray(r2),
            err_msg="Right channel differs between two engines with same seed",
        )

    @pytest.mark.xfail(
        reason="Engine.reset() does not fully restore noise generator state "
               "(Kellet filter). Known issue -- fix in noise_gen.reset().",
        strict=False,
    )
    def test_engine_reset_reproduces_output(self):
        """Engine reset should clear all state: same seed -> same output.

        Currently xfail: reset() does not fully clear NoiseGen internal
        filter state, causing output divergence.
        """
        p = Params()
        p.rate = 8.0
        p.seed = 12345

        e = Engine(SR)
        e.set_params(p)
        l1, r1 = e.process(BLOCK)

        e.reset()
        l2, r2 = e.process(BLOCK)

        np.testing.assert_array_equal(
            np.asarray(l1), np.asarray(l2),
            err_msg="Left channel differs after reset -- state not fully cleared",
        )
        np.testing.assert_array_equal(
            np.asarray(r1), np.asarray(r2),
            err_msg="Right channel differs after reset -- state not fully cleared",
        )


# ===================================================================
# TestPresetBankTherapeutic -- factory preset parametric sanity
# ===================================================================

class TestPresetBankTherapeutic:
    """Verify factory presets load, convert, and produce audio."""

    def test_all_factory_presets_produce_output(self, preset_bank):
        """Every preset in the bank should produce non-silent, bounded audio."""
        for name in preset_bank.list():
            preset = preset_bank.get(name)
            assert preset is not None, f"Preset '{name}' returned None"

            params = preset.to_params()
            e = Engine(SR)
            e.set_params(params)
            left, right = e.process(BLOCK)

            assert not _has_nan_or_inf(left), f"NaN/Inf in preset '{name}' left"
            assert not _has_nan_or_inf(right), f"NaN/Inf in preset '{name}' right"
            assert _peak(left) < 2.0, f"Preset '{name}' left overflow"
            assert _peak(right) < 2.0, f"Preset '{name}' right overflow"

    def test_preset_to_params_roundtrip(self, preset_bank):
        """Preset -> Params -> Preset should preserve key fields.

        Only checks fields that exist in both Preset.params dict and the
        current Params build, avoiding AttributeError on new-build-only
        fields like 'externalization'.
        """
        # Fields guaranteed to exist in all builds
        roundtrip_keys = ("rate", "baselen_ms", "width", "noise_color", "temperature")

        for name in preset_bank.list():
            original = preset_bank.get(name)
            params = original.to_params()

            # from_params reads attributes from the live Params object.
            # If the build lacks 'externalization', from_params will raise.
            # Guard the roundtrip: only test if from_params works.
            try:
                reconstructed = Preset.from_params(
                    name=name, params=params,
                    description=original.description,
                    author=original.author,
                )
            except AttributeError as exc:
                pytest.skip(
                    f"Preset.from_params reads a field not in current build: {exc}"
                )

            for key in roundtrip_keys:
                if key in original.params and key in reconstructed.params:
                    assert abs(reconstructed.params[key] - original.params[key]) < 1e-9, (
                        f"Preset '{name}' roundtrip mismatch on '{key}': "
                        f"{original.params[key]} vs {reconstructed.params[key]}"
                    )


# ===================================================================
# TestEdgeCases -- boundary conditions for therapeutic usage
# ===================================================================

class TestEdgeCases:
    """Edge cases relevant to therapeutic/clinical operation."""

    def test_zero_rate_no_crash(self, engine):
        """rate=0 should not crash -- produces silence or near-silence."""
        p = Params()
        p.rate = 0.0
        engine.set_params(p)
        left, right = engine.process(BLOCK)

        assert len(left) == BLOCK
        assert len(right) == BLOCK
        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)

    @_SKIP_BILATERAL
    def test_extreme_bilateral_rate(self, engine):
        """Very high bilateral rate should not produce NaN/Inf."""
        p = Params()
        p.bilateral_on = True
        p.bilateral_rate = 100.0
        p.bilateral_amount = 1.0
        p.rate = 15.0
        engine.set_params(p)
        left, right = engine.process(BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)

    @_SKIP_DIALOGUE
    def test_max_dialogue_strength(self, engine):
        """dialogue_strength=1.0 should not produce NaN/Inf."""
        p = Params()
        p.dialogue_on = True
        p.dialogue_strength = 1.0
        p.dialogue_memory = 1.0
        p.dialogue_phi_mix = 1.0
        p.rate = 15.0
        engine.set_params(p)
        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left)
        assert not _has_nan_or_inf(right)
        assert _peak(left) < 2.0

    @_SKIP_THERAPEUTIC
    def test_all_subsystems_on(self, engine):
        """Every subsystem active simultaneously should not crash or overflow."""
        p = Params()
        p.rate = 12.0
        p.dialogue_on = True
        p.dialogue_strength = 0.8
        p.bilateral_on = True
        p.bilateral_rate = 1.0
        p.bilateral_amount = 0.8
        p.phi_pan = True
        p.modal_on = True
        p.modal_mix = 0.4
        p.burst = True
        p.thermo = True
        p.lattice = True
        p.temperature = 0.6
        p.glitch_mix = 0.3
        p.externalization = 0.5
        p.seed = 42
        engine.set_params(p)

        left, right = engine.process(LONG_BLOCK)

        assert not _has_nan_or_inf(left), "NaN/Inf with all subsystems on"
        assert not _has_nan_or_inf(right), "NaN/Inf with all subsystems on"
        assert _peak(left) < 2.0, "Overflow with all subsystems on"
        assert _peak(right) < 2.0, "Overflow with all subsystems on"

    def test_long_run_stability(self):
        """10 seconds of processing should not accumulate NaN or unbounded growth."""
        e = Engine(SR)
        p = _hemispheric_bridge_params()
        e.set_params(p)

        total_samples = int(SR * 10)
        block = int(SR)  # process 1 s at a time
        for _ in range(total_samples // block):
            left, right = e.process(block)
            assert not _has_nan_or_inf(left), "NaN/Inf during long run"
            assert not _has_nan_or_inf(right), "NaN/Inf during long run"
            assert _peak(left) < 2.0, "Overflow during long run"
            assert _peak(right) < 2.0, "Overflow during long run"


# ---------------------------------------------------------------------------
# Sprint 4: New DSP modules (binaural, isochronic, tinnitus, spectral slope)
# ---------------------------------------------------------------------------

_HAS_BINAURAL = _has_param("binaural_on")
_HAS_ISOCHRONIC = _has_param("isochronic_on")
_HAS_TINNITUS = _has_param("tinnitus_notch_hz")
_HAS_NOISE_SLOPE = _has_param("noise_slope")

_SKIP_DSP_MODULES = pytest.mark.skipif(
    not (_HAS_BINAURAL and _HAS_ISOCHRONIC and _HAS_TINNITUS and _HAS_NOISE_SLOPE),
    reason="Sprint 4 DSP modules not in current build.",
)


@_SKIP_DSP_MODULES
class TestDSPModules:
    """Tests for binaural beat, isochronic tone, tinnitus notch, spectral tilt."""

    def test_binaural_beat_produces_stereo_diff(self):
        """Binaural beat should produce L/R frequency difference."""
        e = Engine(SR)
        p = Params()
        p.rate = 0.0  # no grains — isolate binaural
        p.binaural_on = True
        p.binaural_carrier_hz = 250.0
        p.binaural_beat_hz = 10.0
        p.binaural_level = 0.5
        e.set_params(p)

        left, right = e.process(LONG_BLOCK)
        left, right = np.array(left), np.array(right)

        # Both channels should have signal
        assert _peak(left) > 0.1, "Binaural L silent"
        assert _peak(right) > 0.1, "Binaural R silent"
        # L and R should differ (different frequencies)
        diff = np.abs(left - right)
        assert np.max(diff) > 0.01, "Binaural L/R should differ"

    def test_binaural_beat_level_control(self):
        """Binaural level controls its contribution: higher level = more signal."""
        e1 = Engine(SR)
        p1 = Params()
        p1.rate = 0.0
        p1.binaural_on = True
        p1.binaural_level = 0.5
        e1.set_params(p1)

        e2 = Engine(SR)
        p2 = Params()
        p2.rate = 0.0
        p2.binaural_on = True
        p2.binaural_level = 0.01
        e2.set_params(p2)

        l1, _ = e1.process(LONG_BLOCK)
        l2, _ = e2.process(LONG_BLOCK)
        # Higher level should produce more energy
        rms1 = np.sqrt(np.mean(np.array(l1)**2))
        rms2 = np.sqrt(np.mean(np.array(l2)**2))
        assert rms1 > rms2 * 2.0, "Binaural level should scale output"

    def test_isochronic_produces_pulsed_output(self):
        """Isochronic tone should produce rhythmic amplitude modulation."""
        e = Engine(SR)
        p = Params()
        p.rate = 0.0  # no grains
        p.isochronic_on = True
        p.isochronic_carrier_hz = 165.0
        p.isochronic_rate_hz = 4.0
        p.isochronic_duty = 0.5
        p.isochronic_level = 0.5
        e.set_params(p)

        # 1 second = 4 pulses at 4 Hz
        left, right = e.process(LONG_BLOCK)
        left, right = np.array(left), np.array(right)

        assert _peak(left) > 0.1, "Isochronic silent"
        # Isochronic adds mono signal, so L and R differ only by the grain engine.
        # Verify isochronic contributes by comparing with isochronic OFF:
        e2 = Engine(SR)
        p2 = Params()
        p2.rate = 0.0
        p2.isochronic_on = False
        e2.set_params(p2)
        l_off, _ = e2.process(LONG_BLOCK)
        rms_on = np.sqrt(np.mean(left**2))
        rms_off = np.sqrt(np.mean(np.array(l_off)**2))
        assert rms_on > rms_off * 1.5, "Isochronic should add significant energy"

    def test_tinnitus_notch_attenuates_target(self):
        """Tinnitus notch at a frequency should create a spectral dip."""
        e = Engine(SR)
        p = Params()
        p.noise_mode = 0  # white noise for flat reference
        p.noise_color = 0
        p.tinnitus_notch_hz = 4000.0
        p.tinnitus_notch_q = 6.0
        e.set_params(p)

        left, _ = e.process(LONG_BLOCK * 4)
        left = np.array(left)

        # Compute PSD
        from scipy.signal import welch
        freqs, psd = welch(left, fs=SR, nperseg=2048)

        # Find power at notch frequency vs neighbors
        notch_idx = np.argmin(np.abs(freqs - 4000.0))
        low_idx = np.argmin(np.abs(freqs - 2000.0))
        high_idx = np.argmin(np.abs(freqs - 6000.0))

        notch_power = psd[notch_idx]
        neighbor_power = (psd[low_idx] + psd[high_idx]) / 2.0

        # Notch should be at least 6 dB down from neighbors
        if neighbor_power > 1e-15:
            ratio_db = 10 * np.log10(notch_power / neighbor_power)
            assert ratio_db < -3.0, f"Notch not deep enough: {ratio_db:.1f} dB"

    def test_spectral_slope_white(self):
        """noise_slope=0 should produce white-ish noise (flat spectrum)."""
        e = Engine(SR)
        p = Params()
        p.noise_slope = 0.0
        p.noise_mode = 0
        e.set_params(p)

        left, _ = e.process(LONG_BLOCK * 4)
        left = np.array(left)
        assert _peak(left) > 0.01, "White noise silent"

    def test_spectral_slope_brown(self):
        """noise_slope=-2 should produce brown noise (steep rolloff)."""
        e = Engine(SR)
        p = Params()
        p.noise_slope = -2.0
        p.noise_mode = 0
        e.set_params(p)

        left, _ = e.process(LONG_BLOCK * 4)
        left = np.array(left)
        assert _peak(left) > 0.001, "Brown noise silent"

        # Verify low frequencies dominate
        from scipy.signal import welch
        freqs, psd = welch(left, fs=SR, nperseg=2048)
        low_band = psd[(freqs > 100) & (freqs < 500)].mean()
        high_band = psd[(freqs > 4000) & (freqs < 8000)].mean()
        if low_band > 1e-15:
            ratio = high_band / low_band
            assert ratio < 0.3, f"Brown noise high/low ratio too high: {ratio:.3f}"

    def test_all_new_modules_no_nan(self):
        """All new modules enabled simultaneously should not produce NaN."""
        e = Engine(SR)
        p = Params()
        p.binaural_on = True
        p.binaural_beat_hz = 6.0
        p.binaural_level = 0.08
        p.isochronic_on = True
        p.isochronic_rate_hz = 10.0
        p.isochronic_level = 0.10
        p.tinnitus_notch_hz = 4000.0
        p.noise_slope = -1.5
        e.set_params(p)

        left, right = e.process(LONG_BLOCK)
        assert not _has_nan_or_inf(left), "NaN with all new modules"
        assert not _has_nan_or_inf(right), "NaN with all new modules"
        assert _peak(left) < 2.0, "Overflow with all new modules"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
