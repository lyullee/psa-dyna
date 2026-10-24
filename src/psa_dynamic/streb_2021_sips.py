"""Pure-component Sips primitives transcribed from Streb & Mazzotti (2021).

This module is deliberately *not* a multicomponent closure.  In particular,
evaluating each species at its partial pressure is not the source's RAST/Wilson
model and must not be used to claim a cyclic replay.  The source contract is
still awaiting complete source-parameter reconciliation and a documented
experimental fit domain.
"""

from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path


SOURCE = Path(__file__).resolve().parent / "literature" / "STREB_MAZZOTTI_2021_SOURCE_CONTRACT_SCREEN.json"
R_KJ_MOL_K = 8.31446261815324e-3


@lru_cache(maxsize=1)
def _source() -> dict:
    data = json.loads(SOURCE.read_text(encoding="utf-8-sig"))
    if data.get("schema") != "dynamic-psa-streb-mazzotti-2021-source-screen/v1":
        raise ValueError("Unexpected Streb source-contract schema")
    if data.get("release_flags", {}).get("engine_input_allowed") is not False:
        raise ValueError("A released engine source must use a separately qualified adapter")
    return data


def parameters(species: str, temperature_k: float) -> tuple[float, float, float]:
    """Return q_inf [mol/kg], b [1/bar], and dimensionless Sips exponent.

    The paper prints T_ref in degrees Celsius, but the ratios and Boltzmann
    factors use absolute kelvin.  B2 and R are both represented in kJ/mol.
    """
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("Temperature must be finite and positive in kelvin")
    data = _source()["pure_component_table_2"]
    if species not in data["components"]:
        raise ValueError(f"Species not in source Table 2: {species}")
    row = data["components"][species]
    tref_k = float(data["reference_temperature_c"]) + 273.15
    qmax = row["A1_mol_kg"] * math.exp(row["A2"] * (temperature_k / tref_k - 1.0))
    affinity = row["B1_bar_inv"] * math.exp(row["B2_kj_mol"] / (R_KJ_MOL_K * temperature_k))
    exponent = row["C1"] + row["C2"] * math.exp(temperature_k / tref_k - 1.0)
    if not all(math.isfinite(value) and value > 0 for value in (qmax, affinity, exponent)):
        raise ValueError("Source Sips parameters are nonfinite or nonpositive at this temperature")
    return qmax, affinity, exponent


def pure_loading_mol_kg(species: str, pressure_bar: float, temperature_k: float) -> float:
    """Evaluate the printed single-component Sips equation only."""
    if not math.isfinite(pressure_bar) or pressure_bar < 0:
        raise ValueError("Pressure must be finite and nonnegative in bar")
    qmax, affinity, exponent = parameters(species, temperature_k)
    if pressure_bar == 0:
        return 0.0
    log_term = exponent * (math.log(affinity) + math.log(pressure_bar))
    if log_term >= 0:
        fraction = 1.0 / (1.0 + math.exp(-log_term))
    else:
        scaled = math.exp(log_term)
        fraction = scaled / (1.0 + scaled)
    return qmax * fraction


def pure_spreading_integral_mol_kg(species: str, pressure_bar: float, temperature_k: float) -> float:
    """Return integral_0^p q(p') d(ln p'), used by the RAST candidate.

    For a Sips isotherm this is q_inf/s * ln(1+(b*p)^s).  This is
    a mathematical pure-component primitive, not a mixture prediction.
    """
    if not math.isfinite(pressure_bar) or pressure_bar < 0:
        raise ValueError("Pressure must be finite and nonnegative in bar")
    qmax, affinity, exponent = parameters(species, temperature_k)
    if pressure_bar == 0:
        return 0.0
    log_term = exponent * (math.log(affinity) + math.log(pressure_bar))
    softplus = max(log_term, 0.0) + math.log1p(math.exp(-abs(log_term)))
    return qmax / exponent * softplus


def inverse_pure_spreading_pressure_bar(species: str, integral_mol_kg: float, temperature_k: float) -> float:
    """Invert the pure Sips spreading integral; zero maps exactly to zero."""
    if not math.isfinite(integral_mol_kg) or integral_mol_kg < 0:
        raise ValueError("Spreading integral must be finite and nonnegative")
    qmax, affinity, exponent = parameters(species, temperature_k)
    if integral_mol_kg == 0:
        return 0.0
    scaled = exponent * integral_mol_kg / qmax
    log_expm1 = (
        math.log(math.expm1(scaled)) if scaled <= 50
        else scaled + math.log1p(-math.exp(-scaled))
    )
    log_pressure = log_expm1 / exponent - math.log(affinity)
    if log_pressure > 700:
        raise ValueError("Inverted pressure exceeds safe finite range")
    return math.exp(log_pressure)
