"""
aureonoise - Phi Head Geometry + Pinna Multipath tests
Verifies the ellipsoid head model, pinna multipath, torso, and distance
response ported from the C++ Max external.
"""

import math
import pytest
from aureonoise._core import PhiModel

PHI = 1.6180339887


def test_geometry_phi_ratio():
    """Vertical axis should be lateral x phi."""
    model = PhiModel(head_b=0.0875, phi_ratio=PHI)
    geom = model.geometry()
    assert abs(geom['axis_vert'] / geom['axis_lat'] - PHI) < 0.001


def test_geometry_nose_extrusion():
    """axis_front = axis_front_base + nose_len."""
    model = PhiModel(head_b=0.0875)
    geom = model.geometry()
    assert abs(geom['axis_front'] - (geom['axis_front_base'] + geom['nose_len'])) < 1e-9


def test_geometry_back_scale():
    """Back axis should be back_scale * head_b."""
    model = PhiModel(head_b=0.0875, back_scale=1.05)
    geom = model.geometry()
    assert abs(geom['axis_back'] - 1.05 * 0.0875) < 1e-9


def test_itd_range():
    """ITD at extreme pan should be within physiological range."""
    model = PhiModel(head_b=0.0875)
    result = model.compute_head(sr=48000.0, pan=1.0, pitch_rad=0.0, phi_distance=0.5)
    itd_us = abs(result['itd_samples']) / 48000.0 * 1e6
    assert 200 < itd_us < 900, f"ITD {itd_us:.0f}us outside physiological range"


def test_itd_zero_at_center():
    """ITD should be ~0 when pan is centered."""
    model = PhiModel(head_b=0.0875)
    result = model.compute_head(sr=48000.0, pan=0.0, pitch_rad=0.0, phi_distance=0.5)
    assert abs(result['itd_samples']) < 1.0, f"ITD at center: {result['itd_samples']}"


def test_itd_symmetry():
    """ITD should be antisymmetric: ITD(pan) = -ITD(-pan)."""
    model = PhiModel(head_b=0.0875)
    left = model.compute_head(sr=48000.0, pan=-0.7, pitch_rad=0.0, phi_distance=0.5)
    right = model.compute_head(sr=48000.0, pan=0.7, pitch_rad=0.0, phi_distance=0.5)
    assert abs(left['itd_samples'] + right['itd_samples']) < 0.01, (
        f"ITD not antisymmetric: {left['itd_samples']} vs {right['itd_samples']}"
    )


def test_itd_distance_scaling():
    """Closer distance should produce slightly larger ITD (itd_scale factor)."""
    model = PhiModel(head_b=0.0875)
    near = model.compute_head(sr=48000.0, pan=1.0, pitch_rad=0.0, phi_distance=0.0)
    far = model.compute_head(sr=48000.0, pan=1.0, pitch_rad=0.0, phi_distance=1.0)
    assert abs(near['itd_samples']) >= abs(far['itd_samples']), (
        f"Near ITD {near['itd_samples']} should >= far ITD {far['itd_samples']}"
    )


def test_pinna_peaks():
    """Pinna should have 3 resonant peaks in expected ranges."""
    model = PhiModel(head_b=0.0875)
    pinna = model.design_pinna(pan=0.5, elev_norm=0.5, elev_notch=0.5)
    # Peak 1 ~ 2900 Hz (ear canal), Peak 2 ~ 4630 Hz (concha), Peak 3 ~ 7800 Hz (cymba)
    assert 2500 < pinna['peak_freqs'][0] < 4000
    assert 3500 < pinna['peak_freqs'][1] < 6000
    assert 6000 < pinna['peak_freqs'][2] < 10000


def test_pinna_tap_count():
    """Pinna should have 5 taps."""
    model = PhiModel(head_b=0.0875)
    pinna = model.design_pinna(pan=0.5, elev_norm=0.0, elev_notch=0.5)
    assert pinna['tap_count'] == 5


def test_pinna_first_tap_zero_delay():
    """First pinna tap (direct path) should have zero delay."""
    model = PhiModel(head_b=0.0875)
    pinna = model.design_pinna(pan=0.5)
    assert pinna['tap_delay_sec'][0] == 0.0


def test_pinna_delays_phi_scaled():
    """Pinna tap delays (k >= 2) should scale by phi."""
    model = PhiModel(head_b=0.0875)
    pinna = model.design_pinna(pan=0.5)
    delays = pinna['tap_delay_sec']
    # taps 2,3,4 should have delta*PHI^n progression
    for k in range(2, 4):
        if delays[k] > 1e-12 and delays[k + 1] > 1e-12:
            ratio = delays[k + 1] / delays[k]
            assert abs(ratio - PHI) < 0.01, (
                f"Delay ratio tap[{k+1}]/tap[{k}] = {ratio:.4f}, expected phi"
            )


def test_pinna_notch_shift():
    """Notch shift should decrease notch frequencies (multiply by 0.90)."""
    model = PhiModel(head_b=0.0875)
    normal = model.design_pinna(pan=0.5, elev_norm=0.5, elev_notch=0.5, notch_shift=False)
    shifted = model.design_pinna(pan=0.5, elev_norm=0.5, elev_notch=0.5, notch_shift=True)
    for i in range(3):
        # Notch freq = c / (2*dd), dd increases with shift => freq decreases
        # Actually: dd *= 0.90 means dd decreases => freq increases
        # Wait: dd *= 0.90 means dd gets smaller, so c/(2*dd) gets larger
        assert shifted['notch_freqs'][i] > normal['notch_freqs'][i], (
            f"Notch shift should increase freq: {shifted['notch_freqs'][i]} vs {normal['notch_freqs'][i]}"
        )


def test_distance_attenuation():
    """Far distance should have lower gain."""
    model = PhiModel(head_b=0.0875)
    near = model.compute_distance(sr=48000.0, phi_distance=0.0)
    far = model.compute_distance(sr=48000.0, phi_distance=1.0)
    assert near['direct_gain'] > far['direct_gain']


def test_distance_hf_rolloff():
    """Far distance should have lower HF weight."""
    model = PhiModel(head_b=0.0875)
    near = model.compute_distance(sr=48000.0, phi_distance=0.0)
    far = model.compute_distance(sr=48000.0, phi_distance=1.0)
    assert near['hf_weight'] > far['hf_weight']


def test_distance_lowpass_alpha():
    """Far distance should have larger lowpass alpha (more filtering)."""
    model = PhiModel(head_b=0.0875)
    near = model.compute_distance(sr=48000.0, phi_distance=0.0)
    far = model.compute_distance(sr=48000.0, phi_distance=1.0)
    assert far['lowpass_alpha'] > near['lowpass_alpha']


def test_torso_inactive_at_zero_mix():
    """Torso should be inactive when mix is 0."""
    model = PhiModel(head_b=0.0875)
    resp = model.compute_torso(sr=48000.0, pan=0.5, torso_mix=0.0)
    assert resp['active'] is False


def test_torso_active_with_mix():
    """Torso should be active when mix > 0."""
    model = PhiModel(head_b=0.0875)
    resp = model.compute_torso(sr=48000.0, pan=0.5, torso_mix=0.3)
    assert resp['active'] is True
    assert resp['tap_count'] == 2


def test_torso_gains_scale_with_mix():
    """Higher torso mix should produce higher gains."""
    model = PhiModel(head_b=0.0875)
    low = model.compute_torso(sr=48000.0, pan=0.5, torso_mix=0.1)
    high = model.compute_torso(sr=48000.0, pan=0.5, torso_mix=0.4)
    assert high['gains'][0] > low['gains'][0]
    assert high['gains'][1] > low['gains'][1]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
