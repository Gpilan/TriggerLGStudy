"""Idealized peak-normalized, two-sided gamma SPE (ns).

The rising gamma branch is time-scaled by rise10-90. The falling branch
has an independent scale fixed by FWHM. Value and first derivative match
at the peak; second derivative need not. This is not a measured PMT fit.
"""
import math
from dataclasses import dataclass
from functools import lru_cache


def _gamma_response(shape, x):
    if x <= 0.0:
        return 0.0
    return math.exp(shape * math.log(x / shape) + shape - x)


def _gamma_root(shape, level, low, high, rising):
    for _ in range(100):
        mid = (low + high) / 2
        if (_gamma_response(shape, mid) < level) == rising:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def _gamma_dimensionless_rise(shape):
    return (_gamma_root(shape, .9, 0., shape, True)
            - _gamma_root(shape, .1, 0., shape, True))


def _gamma_tau_from_rise(response_rise_ns, gamma_shape):
    return response_rise_ns / _gamma_dimensionless_rise(gamma_shape)


def _fall_root(shape, level):
    high = max(2 * shape, 1.)
    while _gamma_response(shape, high) > level:
        high *= 2
    return _gamma_root(shape, level, shape, high, False)


@dataclass(frozen=True)
class Response:
    shape: float
    rise_tau_ns: float
    fall_tau_ns: float
    peak_ns: float
    duration_ns: float
    minimum_fwhm_ns: float

    def value(self, t):
        x = (t / self.rise_tau_ns if t <= self.peak_ns else
             self.shape + (t - self.peak_ns) / self.fall_tau_ns)
        return _gamma_response(self.shape, x)


@lru_cache(maxsize=256)
def resolve_response(rise_ns, fwhm_ns, shape=1.915604733026,
                     tau_ns=None, tail_level=1e-10):
    for name, value in (("rise", rise_ns), ("FWHM", fwhm_ns), ("shape", shape)):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    # Bound the supported numerical family rather than silently losing precision.
    if not 0.1 <= shape <= 100:
        raise ValueError("supported gamma shape range is [0.1, 100]")
    if not math.isfinite(tail_level) or not 0 < tail_level <= 1e-4:
        raise ValueError("kernel tail level must be finite in (0, 1e-4]")
    rise_tau = _gamma_tau_from_rise(rise_ns, shape)
    if tau_ns is not None:
        if not math.isfinite(tau_ns) or not math.isclose(tau_ns, rise_tau, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError("--gamma-tau-ns conflicts with requested rise; omit tau or set a consistent rise")
    peak = shape * rise_tau
    minimum = (shape - _gamma_root(shape, .5, 0., shape, True)) * rise_tau
    if fwhm_ns <= minimum:
        raise ValueError(f"FWHM must exceed {minimum:.9g} ns for rise={rise_ns}, shape={shape}")
    fall_tau = (fwhm_ns - minimum) / (_fall_root(shape, .5) - shape)
    duration = peak + (_fall_root(shape, tail_level) - shape) * fall_tau
    return Response(shape, rise_tau, fall_tau, peak, duration, minimum)


def sample_kernel(response, step_ns):
    if not math.isfinite(step_ns) or step_ns <= 0:
        raise ValueError("sample step must be finite and positive")
    return [response.value(i * step_ns)
            for i in range(int(math.ceil(response.duration_ns / step_ns)) + 1)]


def convolve_impulses(impulse, kernel):
    """Nonnegative linear convolution: sparse shift/add or zero-padded real FFT."""
    import numpy as np
    values=np.asarray(impulse,dtype=float);response=np.asarray(kernel,dtype=float)
    if not len(values) or not len(response):
        raise ValueError("convolution inputs must be nonempty")
    nonzero=np.flatnonzero(values)
    if len(nonzero)<=8:
        voltage=np.zeros(len(values))
        for i in nonzero:
            length=min(len(response),len(values)-i)
            voltage[i:i+length]+=values[i]*response[:length]
    else:
        length=len(values)+len(response)-1
        padded=1 << (length-1).bit_length()
        voltage=np.fft.irfft(np.fft.rfft(values,padded)*np.fft.rfft(response,padded),padded)[:len(values)]
        # For nonnegative photoelectron impulses/SPEs, negative values are roundoff.
        voltage=np.maximum(voltage,0.)
    return voltage.tolist()
