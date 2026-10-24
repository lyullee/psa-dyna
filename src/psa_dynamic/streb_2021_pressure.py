"""Streb--Mazzotti Part 2 pressure-profile mathematics, not an engine driver.

Section 2.4 specifies linear pressurization and exponential BD/BD-vac
boundary pressures.  The published decay coefficients were fitted to the
base cyclic experiment.  This module evaluates those functions only; it does
not impose a Dirichlet pressure on the finite-volume bed or claim that the
valve/MFC diagnostic reproduces the measured trajectory.
"""

from __future__ import annotations

import math


SOURCE_URL = "https://link.springer.com/article/10.1007/s10450-021-00308-w"


def _finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def linear_pressure_bar(start_bar: float, end_bar: float,
                        duration_s: float, elapsed_s: float) -> tuple[float, float]:
    """Return prescribed pressure and dP/dt for an explicitly timed Press.

    The source does not report the split between pressurization and adsorption
    inside Feed.  A caller must supply that duration; this function never
    chooses it or clips time outside the stated interval.
    """
    start = _finite(start_bar, "start pressure")
    end = _finite(end_bar, "end pressure")
    duration = _finite(duration_s, "duration")
    elapsed = _finite(elapsed_s, "elapsed time")
    if min(start, end) <= 0 or duration <= 0 or not 0 <= elapsed <= duration:
        raise ValueError("Linear pressure requires positive endpoints/duration and in-step time")
    slope = (end - start) / duration
    return start + slope * elapsed, slope


def exponential_pressure_bar(start_bar: float, asymptote_bar: float,
                             decay_per_s: float, elapsed_s: float) -> tuple[float, float]:
    """Evaluate Part 2 Eq. (3): P_low+(P_high-P_low) exp(-xi t).

    The second argument is the mathematical asymptote, not automatically the
    measured finite-time endpoint.  Return both P and its exact derivative.
    """
    start = _finite(start_bar, "start pressure")
    low = _finite(asymptote_bar, "asymptotic pressure")
    decay = _finite(decay_per_s, "decay coefficient")
    elapsed = _finite(elapsed_s, "elapsed time")
    if min(start, low) <= 0 or start <= low or decay <= 0 or elapsed < 0:
        raise ValueError("Exponential decay requires start > positive asymptote, xi > 0, t >= 0")
    remainder = (start - low) * math.exp(-decay * elapsed)
    return low + remainder, -decay * remainder


def asymptote_for_finite_endpoint_bar(start_bar: float, end_bar: float,
                                      decay_per_s: float, duration_s: float) -> float:
    """Invert Eq. (3) if a finite-time endpoint is known exactly.

    This derived asymptote is a diagnostic mathematical quantity, *not* a
    published Streb parameter or permission to fit an experimental endpoint.
    """
    start = _finite(start_bar, "start pressure")
    end = _finite(end_bar, "end pressure")
    decay = _finite(decay_per_s, "decay coefficient")
    duration = _finite(duration_s, "duration")
    if start <= end or end <= 0 or decay <= 0 or duration <= 0:
        raise ValueError("Endpoint inversion requires start > end > 0, xi > 0, duration > 0")
    residual = math.exp(-decay * duration)
    asymptote = (end - start * residual) / (1.0 - residual)
    if asymptote <= 0:
        raise ValueError("Requested finite endpoint implies a nonphysical asymptote")
    return asymptote
