"""
aureonoise - Granular noise/glitch texture generator

φ-based stochastic audio synthesis with binaural spatialization.
"""

from aureonoise._core import (
    # Constants
    PHI,
    INV_PHI,
    PHI_SQ,
    INV_PHI_SQ,
    INV_PHI_CU,
    MAX_GRAINS,
    RING_SIZE,
    
    # Core classes
    Rng,
    Weyl,
    NoiseColor,
    NoiseColorState,
    GrainKind,
    RingBuffer,
    OrnsteinUhlenbeck,
    Lattice,
    Hawkes,
    EnvelopeShape,
    Envelope,
    Grain,
    GrainPool,
    Params,
    Engine,
)

__version__ = "0.1.0"
__author__ = "Alessio Ivoy Cazzaniga"
__license__ = "MIT"

__all__ = [
    # Constants
    "PHI",
    "INV_PHI",
    "PHI_SQ",
    "INV_PHI_SQ",
    "INV_PHI_CU",
    "MAX_GRAINS",
    "RING_SIZE",
    
    # Core classes
    "Rng",
    "Weyl",
    "NoiseColor",
    "NoiseColorState",
    "GrainKind",
    "RingBuffer",
    "OrnsteinUhlenbeck",
    "Lattice",
    "Hawkes",
    "EnvelopeShape",
    "Envelope",
    "Grain",
    "GrainPool",
    "Params",
    "Engine",
]
