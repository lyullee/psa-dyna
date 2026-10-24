"""Resolve scalar stream states only from phase-aligned trajectory samples.

The cycle-average PSA ledger contains composition and integrated enthalpy but
does not, by itself, define a thermodynamic pressure or temperature.  This
module accepts an explicit ``phase_aligned_streams`` record produced from
sampled outlet histories.  It never falls back to case boundary values.
"""

from __future__ import annotations

import math
from typing import Any


STREAMS = ("FEED", "PRODUCT", "TAIL")


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _resolve_one(name: str, samples: Any) -> dict[str, Any]:
    if not isinstance(samples, list) or not samples:
        return {"status": "unresolved", "errors": ["samples_missing_or_empty"], "stream": name}
    errors: list[str] = []
    previous_time: float | None = None
    weighted_pressure = 0.0
    weighted_temperature = 0.0
    total_moles = 0.0
    positive_samples = 0
    for index, sample in enumerate(samples):
        if not isinstance(sample, dict):
            errors.append(f"sample_{index}:not_object")
            continue
        time_s = _finite(sample.get("time_s"))
        duration_s = _finite(sample.get("duration_s"))
        flow_mol_s = _finite(sample.get("flow_mol_s"))
        pressure_bar = _finite(sample.get("pressure_bar"))
        temperature_k = _finite(sample.get("temperature_k"))
        if time_s is None or time_s < 0.0:
            errors.append(f"sample_{index}:invalid_time_s")
        elif previous_time is not None and time_s < previous_time:
            errors.append(f"sample_{index}:nonmonotonic_time_s")
        else:
            previous_time = time_s
        if duration_s is None or duration_s <= 0.0:
            errors.append(f"sample_{index}:invalid_duration_s")
        if flow_mol_s is None or flow_mol_s < 0.0:
            errors.append(f"sample_{index}:invalid_flow_mol_s")
        if pressure_bar is None or pressure_bar <= 0.0:
            errors.append(f"sample_{index}:invalid_pressure_bar")
        if temperature_k is None or temperature_k <= 0.0:
            errors.append(f"sample_{index}:invalid_temperature_k")
        if errors and any(item.startswith(f"sample_{index}:") for item in errors):
            continue
        moles = float(flow_mol_s * duration_s)
        if moles <= 0.0:
            continue
        positive_samples += 1
        total_moles += moles
        weighted_pressure += moles * pressure_bar
        weighted_temperature += moles * temperature_k
    if errors:
        return {"status": "unresolved", "errors": sorted(set(errors)), "stream": name}
    if positive_samples == 0 or total_moles <= 0.0:
        return {"status": "unresolved", "errors": ["no_positive_molar_exposure"], "stream": name}
    return {
        "status": "resolved",
        "stream": name,
        "sample_count": len(samples),
        "positive_sample_count": positive_samples,
        "molar_exposure_mol": total_moles,
        "pressure_bar": weighted_pressure / total_moles,
        "temperature_k": weighted_temperature / total_moles,
        "basis": "molar-exposure-weighted phase-aligned samples; no boundary-value fallback",
    }


def resolve_phase_aligned_states(result: dict[str, Any]) -> dict[str, Any]:
    """Resolve FEED/PRODUCT/TAIL scalar states or return explicit blockers."""

    source = result.get("phase_aligned_streams")
    streams = source if isinstance(source, dict) else {}
    resolved = {name: _resolve_one(name, streams.get(name)) for name in STREAMS}
    all_resolved = all(row.get("status") == "resolved" for row in resolved.values())
    return {
        "schema": "dynamic-psa-phase-aligned-state/v1",
        "status": "resolved" if all_resolved else "unresolved",
        "streams": resolved,
        "interpretation": (
            "Scalar states are derived only from explicit phase-aligned stream samples. "
            "Missing or invalid samples keep the state unresolved and block DWSIM transfer."
        ),
    }
