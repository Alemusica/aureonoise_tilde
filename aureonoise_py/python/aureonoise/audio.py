"""
aureonoise - Real-time audio engine wrapper

Provides high-level audio I/O with sounddevice.
"""

import numpy as np
from typing import Callable, Optional
import threading
import queue

try:
    import sounddevice as sd
except ImportError:
    sd = None

from aureonoise import Engine, Params

# Hard safety ceiling — independent of user gain control.
# -1 dBFS = 10^(-1/20) ≈ 0.891. Cannot be bypassed by parameter settings.
MAX_SAFE_AMPLITUDE = 0.891


class AudioEngine:
    """Real-time audio engine with sounddevice backend."""

    def __init__(
        self,
        sample_rate: float = 44100.0,
        block_size: int = 512,
        output_device: Optional[int] = None,
    ):
        """
        Initialize the audio engine.

        Args:
            sample_rate: Sample rate in Hz
            block_size: Audio block size (latency trade-off)
            output_device: Output device index (None = default)
        """
        if sd is None:
            raise ImportError("sounddevice is required: pip install sounddevice")

        self.sample_rate = sample_rate
        self.block_size = block_size
        self.output_device = output_device

        # DSP engine
        self.engine = Engine(sample_rate)
        self._params = Params()

        # Audio state
        self._stream: Optional[sd.OutputStream] = None
        self._running = False
        self._lock = threading.Lock()

        # Metering
        self._peak_l = 0.0
        self._peak_r = 0.0
        self._rms_l = 0.0
        self._rms_r = 0.0

        # Analysis hook (T7.3)
        self._analysis_interval = 0  # blocks between analysis (0 = off)
        self._analysis_counter = 0
        self._analysis_buffer_l: list = []
        self._analysis_buffer_r: list = []
        self._on_analysis: Optional[Callable] = None

        # Callbacks
        self._on_meter: Optional[Callable] = None

        # Session timer — accumulated playback time in seconds
        self._session_frames: int = 0  # total frames rendered (atomic int, safe)
        self._max_session_duration: float = 3600.0  # 1 hour default
        self._session_warning_fired: bool = False
        self._session_limit_fired: bool = False
        self._on_session_warning: Optional[Callable] = None
        self._on_session_limit: Optional[Callable] = None
    
    @property
    def params(self) -> Params:
        """Get current parameters."""
        return self._params
    
    @params.setter
    def params(self, p: Params):
        """Set parameters (thread-safe)."""
        with self._lock:
            self._params = p
            self.engine.set_params(p)
    
    def set_param(self, name: str, value):
        """Set a single parameter by name."""
        with self._lock:
            if hasattr(self._params, name):
                setattr(self._params, name, value)
                self.engine.set_params(self._params)
    
    def get_param(self, name: str):
        """Get a single parameter by name."""
        return getattr(self._params, name, None)
    
    def start(self):
        """Start audio output."""
        if self._running:
            return

        # Use device's native sample rate to avoid PortAudio resampling errors
        actual_sr = self.sample_rate
        try:
            dev_info = sd.query_devices(self.output_device, 'output')
            actual_sr = dev_info['default_samplerate']
        except Exception:
            pass

        if actual_sr != self.sample_rate:
            self.sample_rate = actual_sr
            self.engine = Engine(actual_sr)
            self.engine.set_params(self._params)

        def callback(outdata, frames, time_info, status):
            if status:
                print(f"Audio status: {status}")

            with self._lock:
                left, right = self.engine.process(frames)

            left_np = np.array(left, dtype=np.float32)
            right_np = np.array(right, dtype=np.float32)

            # Hard safety limiter — independent of user gain
            left_np = np.clip(left_np, -MAX_SAFE_AMPLITUDE, MAX_SAFE_AMPLITUDE)
            right_np = np.clip(right_np, -MAX_SAFE_AMPLITUDE, MAX_SAFE_AMPLITUDE)

            # Write to output buffer
            outdata[:, 0] = left_np
            outdata[:, 1] = right_np

            # Session timer — accumulate rendered frames (int add, GIL-safe)
            # Polled from meter thread via session_duration property — no callbacks
            self._session_frames += frames

            # Update meters
            self._peak_l = max(self._peak_l * 0.95, np.abs(left_np).max())
            self._peak_r = max(self._peak_r * 0.95, np.abs(right_np).max())
            self._rms_l = np.sqrt(np.mean(left_np**2))
            self._rms_r = np.sqrt(np.mean(right_np**2))

            if self._on_meter:
                self._on_meter(self._peak_l, self._peak_r, self._rms_l, self._rms_r)

            # Real-time analysis hook (T7.3)
            if self._analysis_interval > 0 and self._on_analysis:
                self._analysis_buffer_l.extend(left_np.tolist())
                self._analysis_buffer_r.extend(right_np.tolist())
                self._analysis_counter += 1
                if self._analysis_counter >= self._analysis_interval:
                    self._analysis_counter = 0
                    try:
                        from aureonoise.analysis import analyze
                        buf_l = np.array(self._analysis_buffer_l)
                        buf_r = np.array(self._analysis_buffer_r)
                        # Collect engine stats for the report
                        eng_stats = None
                        try:
                            eng_stats = {
                                'handshake_rate': self.engine.handshake_ratio(),
                                'coherence': self.engine.coherence_mean(),
                                'plv': self.engine.handshake_plv(),
                            }
                        except Exception:
                            pass
                        report = analyze(buf_l, buf_r, self.sample_rate, eng_stats)
                        self._on_analysis(report)
                    except Exception:
                        pass
                    self._analysis_buffer_l.clear()
                    self._analysis_buffer_r.clear()

        self._stream = sd.OutputStream(
            samplerate=self.sample_rate,
            blocksize=self.block_size,
            device=self.output_device,
            channels=2,
            dtype=np.float32,
            callback=callback,
            latency='high',
        )
        self._stream.start()
        self._running = True
    
    def stop(self):
        """Stop audio output."""
        if self._stream is not None:
            self._stream.stop()
            self._stream.close()
            self._stream = None
        self._running = False
    
    def reset(self):
        """Reset the DSP engine."""
        with self._lock:
            self.engine.reset()
    
    def is_running(self) -> bool:
        """Check if audio is running."""
        return self._running
    
    def on_meter(self, callback: Callable[[float, float, float, float], None]):
        """Set meter callback: fn(peak_l, peak_r, rms_l, rms_r)"""
        self._on_meter = callback
    
    def on_analysis(self, callback: Callable, interval_blocks: int = 86):
        """Set analysis callback: fn(AnalysisReport). Called every interval_blocks.
        Default 86 blocks ≈ 1 second at 512 block size / 44100 Hz."""
        self._on_analysis = callback
        self._analysis_interval = max(1, interval_blocks)

    # ── Session timer ────────────────────────────────────────────────

    @property
    def session_duration(self) -> float:
        """Elapsed playback time in seconds (accumulated from callback frames)."""
        return self._session_frames / self.sample_rate if self.sample_rate > 0 else 0.0

    @property
    def max_session_duration(self) -> float:
        """Maximum session duration in seconds (default 3600 = 1 hour)."""
        return self._max_session_duration

    @max_session_duration.setter
    def max_session_duration(self, value: float):
        self._max_session_duration = max(60.0, value)  # minimum 1 minute

    def on_session_warning(self, callback: Optional[Callable]):
        """Set callback fired at 80% of max session duration. fn(elapsed, max)."""
        self._on_session_warning = callback

    def on_session_limit(self, callback: Optional[Callable]):
        """Set callback fired at 100% of max session duration. fn(elapsed, max).
        Does NOT auto-stop — advisory only."""
        self._on_session_limit = callback

    def reset_session_timer(self):
        """Reset session timer to zero. Call when starting a new session."""
        self._session_frames = 0
        self._session_warning_fired = False
        self._session_limit_fired = False

    def get_meters(self) -> tuple:
        """Get current meter values: (peak_l, peak_r, rms_l, rms_r)"""
        return (self._peak_l, self._peak_r, self._rms_l, self._rms_r)
    
    @staticmethod
    def list_devices():
        """List available audio devices."""
        if sd is None:
            return []
        return sd.query_devices()
    
    @staticmethod
    def default_device():
        """Get default output device index."""
        if sd is None:
            return None
        return sd.default.device[1]  # Output device


def render_offline(
    params: Params,
    duration_sec: float,
    sample_rate: float = 44100.0,
) -> tuple:
    """
    Render audio offline (non-realtime).
    
    Args:
        params: DSP parameters
        duration_sec: Duration in seconds
        sample_rate: Sample rate in Hz
    
    Returns:
        Tuple of (left_channel, right_channel) as numpy arrays
    """
    engine = Engine(sample_rate)
    engine.set_params(params)
    
    num_samples = int(duration_sec * sample_rate)
    left, right = engine.process(num_samples)
    
    return np.array(left), np.array(right)


def save_wav(
    filename: str,
    left: np.ndarray,
    right: np.ndarray,
    sample_rate: int = 44100,
    normalize: bool = True,
):
    """
    Save audio to WAV file.
    
    Args:
        filename: Output filename
        left: Left channel
        right: Right channel  
        sample_rate: Sample rate in Hz
        normalize: Normalize to -0.5 dB peak
    """
    import wave
    import struct
    
    # Interleave
    stereo = np.column_stack([left, right])
    
    if normalize:
        peak = np.abs(stereo).max()
        if peak > 1e-6:
            target = 10**(-0.5/20)  # -0.5 dB
            stereo = stereo * (target / peak)
    
    # Clip
    stereo = np.clip(stereo, -1.0, 1.0)
    
    # Convert to 16-bit
    stereo_int = (stereo * 32767).astype(np.int16)
    
    with wave.open(filename, 'wb') as wav:
        wav.setnchannels(2)
        wav.setsampwidth(2)
        wav.setframerate(sample_rate)
        wav.writeframes(stereo_int.tobytes())
