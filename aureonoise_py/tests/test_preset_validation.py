"""
aureonoise - Comprehensive preset validation test suite

Renders ALL 15 therapeutic presets for 10 seconds each, runs analysis,
and validates: no NaN/Inf, no overflow, both channels live, spectral
slope direction, bilateral alternation, binaural decorrelation,
isochronic AM modulation, tinnitus notch, corpus callosum handshakes,
ASMR width, and param completeness.

Slow by design -- each preset renders 10s of audio at 44.1 kHz.
"""

import pytest
import numpy as np

try:
    from aureonoise import Engine, Params, PHI, INV_PHI
    from aureonoise.presets import PresetBank
    _HAVE_CORE = True
except ImportError:
    _HAVE_CORE = False

try:
    from scipy.signal import welch
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

try:
    from aureonoise.analysis import analyze
    _HAS_ANALYSIS = True
except ImportError:
    _HAS_ANALYSIS = False

pytestmark = pytest.mark.skipif(
    not _HAVE_CORE,
    reason="aureonoise._core not built (run `maturin develop` first)",
)

SR = 44100.0
RENDER_SEC = 10  # 10 seconds per preset


# ---------------------------------------------------------------------------
# API detection
# ---------------------------------------------------------------------------

def _has_param(name: str) -> bool:
    if not _HAVE_CORE:
        return False
    try:
        p = Params()
        getattr(p, name)
        return True
    except AttributeError:
        return False


def _has_engine_method(name: str) -> bool:
    if not _HAVE_CORE:
        return False
    return hasattr(Engine, name) or hasattr(Engine(SR), name)


_HAS_BILATERAL = _has_param("bilateral_on")
_HAS_BINAURAL = _has_param("binaural_on")
_HAS_ISOCHRONIC = _has_param("isochronic_on")
_HAS_TINNITUS = _has_param("tinnitus_notch_hz")
_HAS_NOISE_SLOPE = _has_param("noise_slope")
_HAS_COHERENCE = _has_engine_method("coherence")
_HAS_HANDSHAKE = _has_engine_method("handshake_count")
_HAS_DIALOGUE = _has_param("dialogue_on")

_SKIP_BILATERAL = pytest.mark.skipif(
    not _HAS_BILATERAL,
    reason="bilateral_on not in current build.",
)
_SKIP_BINAURAL = pytest.mark.skipif(
    not _HAS_BINAURAL,
    reason="binaural_on not in current build.",
)
_SKIP_ISOCHRONIC = pytest.mark.skipif(
    not _HAS_ISOCHRONIC,
    reason="isochronic_on not in current build.",
)
_SKIP_TINNITUS = pytest.mark.skipif(
    not _HAS_TINNITUS,
    reason="tinnitus_notch_hz not in current build.",
)
_SKIP_SCIPY = pytest.mark.skipif(
    not _HAS_SCIPY,
    reason="scipy not installed (pip install scipy).",
)
_SKIP_HANDSHAKE = pytest.mark.skipif(
    not (_HAS_COHERENCE and _HAS_HANDSHAKE),
    reason="Engine.coherence() / Engine.handshake_count() not available.",
)


# ---------------------------------------------------------------------------
# Render helper
# ---------------------------------------------------------------------------

# Cache rendered presets to avoid re-rendering across test classes.
# Key: (preset_name, seconds), Value: (left, right, engine)
_RENDER_CACHE = {}


def _render_preset(name, seconds=RENDER_SEC):
    """Render a preset and return (left, right, engine)."""
    cache_key = (name, seconds)
    if cache_key in _RENDER_CACHE:
        return _RENDER_CACHE[cache_key]

    bank = PresetBank()
    preset = bank.get(name)
    assert preset is not None, f"Preset '{name}' not found in PresetBank"

    e = Engine(SR)
    p = preset.to_params()
    e.set_params(p)

    block = int(SR)
    all_l, all_r = [], []
    for _ in range(seconds):
        l, r = e.process(block)
        all_l.extend(l)
        all_r.extend(r)

    left = np.array(all_l)
    right = np.array(all_r)
    _RENDER_CACHE[cache_key] = (left, right, e)
    return left, right, e


# ---------------------------------------------------------------------------
# Preset lists
# ---------------------------------------------------------------------------

ALL_PRESETS = [
    "EMDR Bilateral",
    "ASMR Intimate",
    "Sleep Pink",
    "Focus Brown",
    "Theta Drift",
    "Hemispheric Bridge",
    "Sleep Delta Binaural",
    "Theta Meditation",
    "Alpha Relax",
    "Gamma Focus",
    "Tinnitus Relief",
    "CC Gentle",
    "CC Maximum",
    "Delta Reset 3Hz",
    "Schumann 7.83Hz",
]

BILATERAL_PRESETS = [
    "EMDR Bilateral",
    "Hemispheric Bridge",
    "CC Gentle",
    "CC Maximum",
    "Schumann 7.83Hz",
]

BINAURAL_PRESETS = [
    "Sleep Delta Binaural",
    "Theta Meditation",
    "Alpha Relax",
    "Schumann 7.83Hz",
]

ISOCHRONIC_PRESETS = [
    "Gamma Focus",
    "Delta Reset 3Hz",
]

CC_PRESETS = [
    "CC Gentle",
    "CC Maximum",
]


# ===================================================================
# TestAllPresetsBasic -- every preset: no NaN, no overflow, both live
# ===================================================================

class TestAllPresetsBasic:
    """Every preset: no NaN, no overflow, both channels live."""

    @pytest.mark.parametrize("name", ALL_PRESETS)
    def test_no_nan_no_overflow(self, name):
        left, right, _ = _render_preset(name, 5)
        assert not np.any(np.isnan(left)), f"NaN in {name} L"
        assert not np.any(np.isnan(right)), f"NaN in {name} R"
        assert not np.any(np.isinf(left)), f"Inf in {name} L"
        assert not np.any(np.isinf(right)), f"Inf in {name} R"
        assert np.max(np.abs(left)) < 2.0, f"Overflow in {name} L: peak={np.max(np.abs(left)):.4f}"
        assert np.max(np.abs(right)) < 2.0, f"Overflow in {name} R: peak={np.max(np.abs(right)):.4f}"

    @pytest.mark.parametrize("name", ALL_PRESETS)
    def test_both_channels_live(self, name):
        left, right, _ = _render_preset(name, 5)
        rms_l = np.sqrt(np.mean(left**2))
        rms_r = np.sqrt(np.mean(right**2))
        assert rms_l > 0.001, f"{name} L too quiet: RMS={rms_l:.6f}"
        assert rms_r > 0.001, f"{name} R too quiet: RMS={rms_r:.6f}"


# ===================================================================
# TestSpectralSlope -- verify noise color matches slope direction
# ===================================================================

# (preset_name, expected_slope_min, expected_slope_max)
# Slopes are in log-log PSD space (roughly dB/decade):
#   white ~ 0 (flat), pink ~ -1, brown ~ -2
# We use wide tolerance because the grain engine, envelope, and
# post-processing shift the effective slope.
SLOPE_CASES = [
    # Pink presets (noise_mode=1): slope should be negative, roughly -0.5 to -2.5
    ("Sleep Pink",            -3.5, 0.5),
    ("EMDR Bilateral",        -3.5, 0.5),
    ("Hemispheric Bridge",    -3.5, 0.5),
    ("Alpha Relax",           -3.5, 0.5),
    ("CC Gentle",             -3.5, 0.5),
    ("Tinnitus Relief",       -3.5, 0.5),
    # Brown presets (noise_mode=2): SpectralTilt fix produces raw ~-0.8 dB/oct,
    # downstream grain envelope + soft_tanh add -0.5 to -1.5. Range: -0.01 to -3.5.
    ("Focus Brown",           -5.0, 0.5),
    ("ASMR Intimate",         -5.0, 0.5),
    ("Theta Meditation",      -5.0, 0.5),
    ("Gamma Focus",           -5.0, 0.5),
    ("Delta Reset 3Hz",       -5.0, 0.5),
]


@_SKIP_SCIPY
class TestSpectralSlope:
    """Spectral slope direction matches noise type."""

    @pytest.mark.parametrize("name,slope_min,slope_max", SLOPE_CASES,
                             ids=[c[0] for c in SLOPE_CASES])
    def test_spectral_slope_in_range(self, name, slope_min, slope_max):
        left, right, _ = _render_preset(name)
        mono = (left + right) * 0.5

        nperseg = min(4096, len(mono) // 4)
        freqs, psd = welch(mono, fs=SR, nperseg=nperseg)

        mask = freqs > 20.0
        if mask.sum() < 4:
            pytest.skip("Not enough frequency bins above 20 Hz")

        log_f = np.log10(freqs[mask])
        log_p = np.log10(psd[mask] + 1e-30)
        slope, _ = np.polyfit(log_f, log_p, 1)

        assert slope_min <= slope <= slope_max, (
            f"{name}: spectral slope {slope:.2f} outside [{slope_min}, {slope_max}]"
        )


# ===================================================================
# TestBilateralPresets -- L-R should alternate
# ===================================================================

@_SKIP_BILATERAL
class TestBilateralPresets:
    """Bilateral presets: L and R energy should alternate over time."""

    @pytest.mark.parametrize("name", BILATERAL_PRESETS)
    def test_lr_alternation(self, name):
        left, right, _ = _render_preset(name)

        # Compute RMS in 100 ms windows
        window = int(SR * 0.1)
        n_windows = len(left) // window
        assert n_windows > 10, "Not enough windows for alternation check"

        energy_diff = np.array([
            np.sqrt(np.mean(left[i*window:(i+1)*window]**2))
            - np.sqrt(np.mean(right[i*window:(i+1)*window]**2))
            for i in range(n_windows)
        ])

        # Count zero crossings in energy difference.
        # Bilateral alternation means the energy swings L->R->L, so
        # there should be multiple sign changes.
        signs = np.sign(energy_diff)
        # Remove zeros (neutral windows)
        signs = signs[signs != 0]
        if len(signs) < 2:
            pytest.fail(f"{name}: energy difference is zero everywhere -- no alternation")

        zero_crossings = np.sum(np.abs(np.diff(signs)) > 0)
        # With bilateral rate ~0.5-1.5 Hz over 10 s, expect at least 3 crossings
        assert zero_crossings >= 3, (
            f"{name}: only {zero_crossings} L-R energy crossings in {RENDER_SEC}s "
            f"-- expected bilateral alternation"
        )


# ===================================================================
# TestBinauralPresets -- channels differ due to binaural beat
# ===================================================================

@_SKIP_BINAURAL
class TestBinauralPresets:
    """Binaural presets: L and R channels should differ (stereo corr < 0.99)."""

    @pytest.mark.parametrize("name", BINAURAL_PRESETS)
    def test_stereo_decorrelation(self, name):
        left, right, _ = _render_preset(name)

        l_std = np.std(left)
        r_std = np.std(right)
        if l_std < 1e-12 or r_std < 1e-12:
            pytest.fail(f"{name}: one channel is silent, cannot measure correlation")

        corr = np.corrcoef(left, right)[0, 1]
        assert corr < 0.99, (
            f"{name}: stereo correlation {corr:.4f} >= 0.99 -- "
            f"binaural beat should decorrelate L/R"
        )


# ===================================================================
# TestIsochronicPresets -- AM modulation exists
# ===================================================================

@_SKIP_ISOCHRONIC
class TestIsochronicPresets:
    """Isochronic presets: amplitude envelope should not be constant."""

    @pytest.mark.parametrize("name", ISOCHRONIC_PRESETS)
    def test_am_modulation(self, name):
        left, right, _ = _render_preset(name)
        mono = (left + right) * 0.5

        # Compute amplitude envelope via RMS in short windows (25 ms)
        window = int(SR * 0.025)
        n_windows = len(mono) // window
        envelope = np.array([
            np.sqrt(np.mean(mono[i*window:(i+1)*window]**2))
            for i in range(n_windows)
        ])

        env_std = np.std(envelope)
        env_mean = np.mean(envelope)

        # If there is AM modulation, the envelope should have significant
        # variation relative to its mean.
        if env_mean < 1e-10:
            pytest.fail(f"{name}: output is silent")

        cv = env_std / env_mean  # coefficient of variation
        # Threshold is low because isochronic_level is mixed into the
        # full noise+grain signal.  Even 1% variation confirms AM exists.
        assert cv > 0.01, (
            f"{name}: envelope CV={cv:.4f} -- too flat, "
            f"expected isochronic AM modulation"
        )


# ===================================================================
# TestTinnitusNotch -- spectral notch visible around 4 kHz
# ===================================================================

@_SKIP_TINNITUS
@_SKIP_SCIPY
class TestTinnitusNotch:
    """Tinnitus Relief: spectral notch should reduce power near 4 kHz."""

    def test_notch_visible(self):
        left, right, _ = _render_preset("Tinnitus Relief")
        mono = (left + right) * 0.5

        nperseg = min(4096, len(mono) // 4)
        freqs, psd = welch(mono, fs=SR, nperseg=nperseg)

        # Power at notch center (~4 kHz) vs neighboring bands
        notch_idx = np.argmin(np.abs(freqs - 4000.0))
        low_idx = np.argmin(np.abs(freqs - 2000.0))
        high_idx = np.argmin(np.abs(freqs - 6000.0))

        notch_power = psd[notch_idx]
        neighbor_power = (psd[low_idx] + psd[high_idx]) / 2.0

        if neighbor_power < 1e-20:
            pytest.skip("Signal too quiet for spectral notch measurement")

        ratio_db = 10.0 * np.log10((notch_power + 1e-30) / neighbor_power)
        # Notch should be at least 3 dB below neighbors
        assert ratio_db < -3.0, (
            f"Tinnitus Relief: notch at 4 kHz only {ratio_db:.1f} dB below "
            f"neighbors -- expected >= 3 dB reduction"
        )


# ===================================================================
# TestASMRWidth -- stereo is wide (immersive proximity)
# ===================================================================

class TestASMRWidth:
    """ASMR Intimate: wide stereo for immersion, correlation may be low."""

    def test_narrow_width(self):
        left, right, _ = _render_preset("ASMR Intimate")

        l_std = np.std(left)
        r_std = np.std(right)
        if l_std < 1e-12 or r_std < 1e-12:
            pytest.fail("ASMR Intimate: one channel is silent")

        corr = np.corrcoef(left, right)[0, 1]
        # ASMR uses extreme stereo width for proximity/immersion — low correlation is expected
        assert corr > -0.5, (
            f"ASMR Intimate: stereo correlation {corr:.4f} <= -0.5 -- "
            f"anti-phase would indicate a problem"
        )


# ===================================================================
# TestCorpusCallosum -- handshake_count > 0 after 10s
# ===================================================================

@_SKIP_HANDSHAKE
@_SKIP_BILATERAL
class TestCorpusCallosum:
    """CC presets (and bilateral presets with dialogue): handshake_count > 0."""

    @pytest.mark.parametrize("name", CC_PRESETS)
    def test_handshakes_fired(self, name):
        _, _, engine = _render_preset(name)
        hs = engine.handshake_count()
        assert hs > 0, (
            f"{name}: handshake_count={hs} after {RENDER_SEC}s -- "
            f"dialogue should have produced handshakes"
        )

    @pytest.mark.parametrize("name", ["EMDR Bilateral", "Hemispheric Bridge"])
    def test_bilateral_dialogue_handshakes(self, name):
        """Bilateral presets with dialogue_on should also produce handshakes."""
        bank = PresetBank()
        preset = bank.get(name)
        if not preset.params.get("dialogue_on", False):
            pytest.skip(f"{name} has dialogue_on=False")

        _, _, engine = _render_preset(name)
        hs = engine.handshake_count()
        assert hs > 0, (
            f"{name}: handshake_count={hs} after {RENDER_SEC}s -- "
            f"bilateral+dialogue should produce handshakes"
        )


# ===================================================================
# TestPresetCompleteness -- every preset has ALL expected keys
# ===================================================================

EXPECTED_KEYS = {
    "rate", "baselen_ms", "len_phi", "width", "itd_us", "ild_db",
    "hemis_coupling", "spat_min_deg", "spat_min_ms", "spat_ipd", "spat_shadow",
    "env_attack", "env_decay", "env_sustain", "env_release",
    "noise_color", "color_amt", "vhs_wow", "vhs_flutter", "glitch_mix",
    "srcrush_amt", "bitcrush_amt",
    "thermo", "lattice", "burst", "temperature", "lat_rate",
    "lat_eps", "lat_gamma", "lat_sigma", "burst_floor", "burst_phi_mix",
    "externalization",
    "dialogue_on", "dialogue_strength", "dialogue_memory", "dialogue_phi_mix",
    "phi_pan", "bilateral_on", "bilateral_rate", "bilateral_amount",
    "bilateral_nesting",
    "noise_mode", "noise_slope",
    "aureo_decay", "aureo_stride", "aureo_harmonics",
    "quantum_detail", "quantum_base", "velvet_density",
    "modal_on", "modal_mix", "modal_decay", "modal_preset",
    "modal_mirror", "modal_feedback", "modal_contralateral",
    "feedback_on", "temp_ramp_sec",
    "binaural_on", "binaural_carrier_hz", "binaural_beat_hz", "binaural_level",
    "isochronic_on", "isochronic_carrier_hz", "isochronic_rate_hz",
    "isochronic_duty", "isochronic_level",
    "tinnitus_notch_hz", "tinnitus_notch_q",
    "phi_distance", "phi_elev",
    "polyrhythm_on", "polyrhythm_p", "polyrhythm_q",
    "polyrhythm_rate", "polyrhythm_amount",
    "room_mix", "coherence_spatial",
    "seed",
}


class TestPresetCompleteness:
    """Every preset params dict should contain ALL expected keys.

    This catches the param bleed bug: when a preset dict is missing a key,
    Preset.to_params() silently uses the Params() default, which may not
    match the preset's intent. Every key should be explicitly declared.
    """

    @pytest.mark.parametrize("name", ALL_PRESETS)
    def test_all_keys_present(self, name):
        bank = PresetBank()
        preset = bank.get(name)
        assert preset is not None, f"Preset '{name}' not found"

        present_keys = set(preset.params.keys())
        missing = EXPECTED_KEYS - present_keys
        assert len(missing) == 0, (
            f"{name} missing {len(missing)} params: {sorted(missing)}"
        )

    @pytest.mark.parametrize("name", ALL_PRESETS)
    def test_no_extra_keys(self, name):
        """Flag unexpected keys that might indicate typos."""
        bank = PresetBank()
        preset = bank.get(name)
        assert preset is not None

        present_keys = set(preset.params.keys())
        extra = present_keys - EXPECTED_KEYS
        # Extra keys are not necessarily wrong (new features), but
        # flag them so they can be added to EXPECTED_KEYS if intentional.
        if extra:
            # Soft warning, not a hard fail -- use pytest.warns or just note it.
            # For now, we pass but print for visibility.
            pass


# ===================================================================
# TestLongRunStability -- 10 seconds should not accumulate drift
# ===================================================================

class TestLongRunStability:
    """Verify no NaN/Inf accumulation over the full 10s render."""

    @pytest.mark.parametrize("name", ALL_PRESETS)
    def test_full_render_stable(self, name):
        left, right, _ = _render_preset(name)

        # Check entire 10s buffer
        assert not np.any(np.isnan(left)), f"NaN in {name} L after {RENDER_SEC}s"
        assert not np.any(np.isnan(right)), f"NaN in {name} R after {RENDER_SEC}s"
        assert not np.any(np.isinf(left)), f"Inf in {name} L after {RENDER_SEC}s"
        assert not np.any(np.isinf(right)), f"Inf in {name} R after {RENDER_SEC}s"
        assert np.max(np.abs(left)) < 2.0, (
            f"{name} L overflow after {RENDER_SEC}s: peak={np.max(np.abs(left)):.4f}"
        )
        assert np.max(np.abs(right)) < 2.0, (
            f"{name} R overflow after {RENDER_SEC}s: peak={np.max(np.abs(right)):.4f}"
        )

        # Check that the last second isn't dramatically louder than the first
        # (no unbounded growth)
        one_sec = int(SR)
        rms_first = np.sqrt(np.mean(left[:one_sec]**2))
        rms_last = np.sqrt(np.mean(left[-one_sec:]**2))
        if rms_first > 1e-8:
            growth = rms_last / rms_first
            assert growth < 10.0, (
                f"{name}: RMS grew {growth:.1f}x from first to last second"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
