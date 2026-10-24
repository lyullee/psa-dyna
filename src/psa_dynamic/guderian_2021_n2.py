"""Source-local constitutive laws for Marcinek et al. (2021) N2 PSA.

The module contains only the O2/N2/CT-350 equations and parameters printed in
doi:10.1007/s10450-021-00320-0.  It deliberately does not reuse H2-PSA or
biogas parameters.  Pressures are in bar, temperatures in kelvin, loadings in
mol/kg, and rates in mol/(kg s) unless a function says otherwise.

The article reports pure temperature-dependent Sips isotherms, classical IAST,
Darken-corrected Arrhenius diffusivity and a QDF rate law.  The public source
does not provide a numeric CSS norm or synchronized raw cycle histories; this
module therefore implements constitutive mathematics, not a validation claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import brentq


R = 8.314462618
REFERENCE_TEMPERATURE_K = 293.15
# Pure isotherms were measured at 293.15--333.15 K, while the same published
# dynamic model was executed against PSA experiments down to 5 C.  Permit
# equation evaluation over that documented model-execution interval, but keep
# the narrower measured equilibrium-fit interval in each case's independent
# trajectory qualification gate.
MIN_TEMPERATURE_K = 278.15
MAX_TEMPERATURE_K = 333.15
# The article states 0--20 bar; the public raw tables contain instrument
# endpoints up to 20.035 bar.  The 0.05-bar allowance includes those measured
# rows without opening a materially wider extrapolation domain.
MAX_PURE_PRESSURE_BAR = 20.05
PARTICLE_RADIUS_M = 0.830e-3
SPECIES_ORDER = ("O2", "N2")

# Chapman--Enskog Lennard--Jones values used for the source Eq. 3 transport
# adapter.  Diameters and epsilon/k are the standard Poling/Bird values also
# tabulated by the open references recorded in the source audit.  They are
# isolated here and are not transferred to H2 or CO2/CH4 branches.
_LJ_SIGMA_ANGSTROM = {"O2": 3.467, "N2": 3.798}
_LJ_EPSILON_OVER_K_K = {"O2": 106.7, "N2": 71.4}
_MOLAR_MASS_G_MOL = {"O2": 31.9988, "N2": 28.0134}


class GuderianSourceDomainError(ValueError):
    """Raised when a requested state exceeds the published source domain."""


@dataclass(frozen=True)
class TemperatureSips:
    saturation_mol_kg: float
    affinity_per_bar: float
    pressure_exponent: float  # article Eq. 8 uses 1/n


@dataclass(frozen=True)
class SourceComponent:
    name: str
    saturation_reference_mol_kg: float
    saturation_temperature_coefficient: float
    affinity_reference_per_bar: float
    affinity_energy_j_mol: float
    n_reference: float
    n_temperature_coefficient: float
    diffusion_frequency_m2_s: float
    activation_energy_j_mol: float
    heat_polynomial: tuple[float, float, float, float, float, float]
    heat_temperature_exponents_k: tuple[float, float, float, float, float, float]


# Table 1 reports kmol/kg.  Multiplication by 1000 gives mol/kg.
# kJ/kmol is numerically identical to J/mol for Q and Ea.
COMPONENTS = {
    "O2": SourceComponent(
        "O2", 3.384, 1.104, 9.436e-2, 1.222e4, 1.120, 3.341e-1,
        1.208e-6, 1622.33,
        (2.384e14, -1.474e12, 4.914e9, -8.498e6, 8.454e3, -1.490e1),
        (-1.084e1, 6.320e1, 4.032e1, 2.323, -2.619e1, 2.584e1),
    ),
    "N2": SourceComponent(
        "N2", 2.707, 1.146, 1.205e-1, 1.187e4, 1.185, 2.263e-1,
        8.186e-5, 25098.55,
        (3.906e14, -1.826e12, 5.312e9, -7.997e6, 6.931e3, -1.363e1),
        (-1.087e1, 8.365e1, 4.669e1, 2.793, -1.921e1, 2.548e1),
    ),
}


def _component(name: str) -> SourceComponent:
    try:
        return COMPONENTS[name]
    except KeyError as error:
        raise ValueError(f"Unsupported Guderian component: {name}") from error


def _temperature_guard(temperature_k: float) -> float:
    temperature = float(temperature_k)
    if not MIN_TEMPERATURE_K <= temperature <= MAX_TEMPERATURE_K:
        raise GuderianSourceDomainError(
            f"Guderian equations are restricted to the published {MIN_TEMPERATURE_K}--{MAX_TEMPERATURE_K} K execution interval"
        )
    return temperature


def temperature_sips(component: str, temperature_k: float) -> TemperatureSips:
    """Evaluate article Eqs. 9--11 using the 20 C reference table."""
    source = _component(component)
    temperature = _temperature_guard(temperature_k)
    saturation = source.saturation_reference_mol_kg * np.exp(
        source.saturation_temperature_coefficient
        * (1.0 - temperature / REFERENCE_TEMPERATURE_K)
    )
    affinity = source.affinity_reference_per_bar * np.exp(
        source.affinity_energy_j_mol / (R * REFERENCE_TEMPERATURE_K)
        * (REFERENCE_TEMPERATURE_K / temperature - 1.0)
    )
    exponent = (
        1.0 / source.n_reference
        + source.n_temperature_coefficient
        * (1.0 - REFERENCE_TEMPERATURE_K / temperature)
    )
    if saturation <= 0.0 or affinity <= 0.0 or exponent <= 0.0:
        raise ArithmeticError("Guderian Sips equation produced non-positive parameters")
    return TemperatureSips(float(saturation), float(affinity), float(exponent))


def pure_loading_mol_kg(component: str, pressure_bar: float, temperature_k: float) -> float:
    """Return pure loading from article Eq. 8 within the 0--20 bar fit range."""
    pressure = float(pressure_bar)
    if pressure < 0.0 or pressure > MAX_PURE_PRESSURE_BAR:
        raise GuderianSourceDomainError("Pure Guderian Sips pressure must be within 0--20.05 bar")
    if pressure == 0.0:
        return 0.0
    sips = temperature_sips(component, temperature_k)
    z = (sips.affinity_per_bar * pressure) ** sips.pressure_exponent
    return float(sips.saturation_mol_kg * z / (1.0 + z))


def spreading_pressure(component: str, pressure_bar: float, temperature_k: float) -> float:
    """Analytic integral of q(P)/P for the source Sips convention."""
    pressure = float(pressure_bar)
    if pressure < 0.0 or pressure > MAX_PURE_PRESSURE_BAR:
        raise GuderianSourceDomainError("IAST pure pressure must be within 0--20.05 bar")
    if pressure == 0.0:
        return 0.0
    sips = temperature_sips(component, temperature_k)
    z = (sips.affinity_per_bar * pressure) ** sips.pressure_exponent
    return float(sips.saturation_mol_kg / sips.pressure_exponent * np.log1p(z))


def pressure_at_spreading_pressure(component: str, value: float, temperature_k: float) -> float:
    """Invert the source Sips spreading-pressure integral analytically."""
    if value < 0.0:
        raise ValueError("Spreading pressure must be non-negative")
    if value == 0.0:
        return 0.0
    sips = temperature_sips(component, temperature_k)
    z = np.expm1(float(value) * sips.pressure_exponent / sips.saturation_mol_kg)
    pressure = z ** (1.0 / sips.pressure_exponent) / sips.affinity_per_bar
    if pressure > MAX_PURE_PRESSURE_BAR * (1.0 + 1e-10):
        raise GuderianSourceDomainError("IAST hypothetical pure pressure exceeds 20.05 bar")
    return float(min(pressure, MAX_PURE_PRESSURE_BAR))


def iast_loading_mol_kg(
    partial_pressures_bar: Sequence[float],
    temperature_k: float,
    species_order: Sequence[str] = SPECIES_ORDER,
) -> np.ndarray:
    """Return source-local classical-IAST component loadings.

    The source reports 11.2, 1.84 and 5.07 percent maximum mixture-loading
    errors at 20, 40 and 60 C in the practical 1--10 bar range.  The result is
    therefore the article model, not exact mixture thermodynamic truth.
    """
    temperature = _temperature_guard(temperature_k)
    names = tuple(species_order)
    if names != SPECIES_ORDER:
        raise ValueError(f"Guderian source order must be {SPECIES_ORDER}")
    partial = np.asarray(partial_pressures_bar, dtype=float)
    if partial.shape != (2,) or np.any(~np.isfinite(partial)) or np.any(partial < 0.0):
        raise ValueError("O2/N2 partial pressures must be two finite non-negative values")
    total = float(partial.sum())
    if total <= 0.0:
        return np.zeros(2)
    if total > MAX_PURE_PRESSURE_BAR:
        raise GuderianSourceDomainError("Guderian IAST total pressure exceeds 20.05 bar")
    active = np.flatnonzero(partial > 0.0)
    if len(active) == 1:
        result = np.zeros(2)
        result[active[0]] = pure_loading_mol_kg(names[active[0]], total, temperature)
        return result
    gas_fraction = partial / total
    maximum_common = min(
        spreading_pressure(names[index], MAX_PURE_PRESSURE_BAR, temperature)
        for index in active
    )

    def allocation_sum(common_spreading_pressure: float) -> float:
        return float(sum(
            gas_fraction[index] * total
            / pressure_at_spreading_pressure(names[index], common_spreading_pressure, temperature)
            for index in active
        ))

    lower = max(1e-14, maximum_common * 1e-12)
    if allocation_sum(maximum_common) > 1.0:
        raise GuderianSourceDomainError(
            "IAST allocation requires a hypothetical pure pressure above 20.05 bar"
        )
    common = brentq(lambda value: allocation_sum(value) - 1.0, lower, maximum_common, xtol=1e-12)
    pure_pressures = np.array([
        pressure_at_spreading_pressure(name, common, temperature) for name in names
    ])
    adsorbed_fraction = gas_fraction * total / pure_pressures
    adsorbed_fraction /= adsorbed_fraction.sum()
    pure_loadings = np.array([
        pure_loading_mol_kg(name, pressure, temperature)
        for name, pressure in zip(names, pure_pressures)
    ])
    total_loading = 1.0 / np.sum(adsorbed_fraction / pure_loadings)
    return adsorbed_fraction * total_loading


def heat_of_adsorption_magnitude_j_mol(
    component: str, loading_mol_kg: float, temperature_k: float
) -> float:
    """Return the positive magnitude plotted in article Figure 2.

    Table 2's polynomial evaluates the signed adsorption enthalpy.  The paper
    plots ``-Delta H`` as a positive magnitude, hence the explicit minus sign.
    The polynomial loading variable is kmol/kg and MJ/kmol converts to J/mol
    by multiplication by 1000.
    """
    source = _component(component)
    temperature = _temperature_guard(temperature_k)
    loading = float(loading_mol_kg)
    sips = temperature_sips(component, temperature)
    if loading < 0.0 or loading > sips.saturation_mol_kg * (1.0 + 1e-10):
        raise GuderianSourceDomainError("Heat correlation loading is outside the Sips capacity")
    w = loading / 1000.0
    powers = (5, 4, 3, 2, 1, 0)
    signed_mj_kmol = sum(
        coefficient * w**power * np.exp(exponent / temperature)
        for coefficient, exponent, power in zip(
            source.heat_polynomial, source.heat_temperature_exponents_k, powers
        )
    )
    magnitude = -float(signed_mj_kmol) * 1000.0
    if not 5.0e3 <= magnitude <= 20.0e3:
        raise GuderianSourceDomainError("Heat correlation left the plotted 5--20 kJ/mol range")
    return magnitude


def adsorption_energy_magnitude_j_kg(
    component: str, loading_mol_kg: float, temperature_k: float
) -> float:
    """Integrate the Table-2 heat magnitude from zero to ``loading``.

    This analytic primitive lets the common engine use the same constitutive
    law in its stored-energy ledger and its differential heat-release term.
    """
    source = _component(component)
    temperature = _temperature_guard(temperature_k)
    loading = float(loading_mol_kg)
    sips = temperature_sips(component, temperature)
    if loading < 0.0 or loading > sips.saturation_mol_kg * (1.0 + 1e-10):
        raise GuderianSourceDomainError("Adsorption-energy loading is outside the Sips capacity")
    w = loading / 1000.0
    powers = (5, 4, 3, 2, 1, 0)
    energy = -1.0e6 * sum(
        coefficient * np.exp(exponent / temperature) * w ** (power + 1) / (power + 1)
        for coefficient, exponent, power in zip(
            source.heat_polynomial, source.heat_temperature_exponents_k, powers
        )
    )
    if not np.isfinite(energy) or energy < 0.0:
        raise GuderianSourceDomainError("Integrated adsorption-energy magnitude is non-positive")
    return float(energy)


def adsorption_energy_temperature_derivative_j_kg_k(
    component: str, loading_mol_kg: float, temperature_k: float
) -> float:
    """Partial derivative of the integrated heat magnitude at fixed loading."""
    source = _component(component)
    temperature = _temperature_guard(temperature_k)
    loading = float(loading_mol_kg)
    sips = temperature_sips(component, temperature)
    if loading < 0.0 or loading > sips.saturation_mol_kg * (1.0 + 1e-10):
        raise GuderianSourceDomainError("Adsorption-energy derivative is outside the Sips capacity")
    w = loading / 1000.0
    powers = (5, 4, 3, 2, 1, 0)
    derivative = 1.0e6 * sum(
        coefficient * exponent * np.exp(exponent / temperature)
        * w ** (power + 1) / ((power + 1) * temperature**2)
        for coefficient, exponent, power in zip(
            source.heat_polynomial, source.heat_temperature_exponents_k, powers
        )
    )
    if not np.isfinite(derivative):
        raise GuderianSourceDomainError("Adsorption-energy temperature derivative is non-finite")
    return float(derivative)


def qdf_rate_mol_kg_s(
    current_loading_mol_kg: Sequence[float],
    equilibrium_loading_mol_kg: Sequence[float],
    temperature_k: float,
) -> np.ndarray:
    """Evaluate article Eqs. 4--5 for O2 and N2 in source order."""
    temperature = _temperature_guard(temperature_k)
    current = np.asarray(current_loading_mol_kg, dtype=float)
    target = np.asarray(equilibrium_loading_mol_kg, dtype=float)
    if current.shape != (2,) or target.shape != (2,) or np.any(~np.isfinite(current)) or np.any(~np.isfinite(target)):
        raise ValueError("QDF loadings must be two finite O2/N2 values")
    if np.any(current <= 0.0) or np.any(target < 0.0):
        raise GuderianSourceDomainError("Published QDF form requires current loading > 0 and target loading >= 0")
    rate = np.empty(2)
    for index, name in enumerate(SPECIES_ORDER):
        source = _component(name)
        sips = temperature_sips(name, temperature)
        theta = current[index] / sips.saturation_mol_kg
        if not 0.0 < theta < 1.0:
            raise GuderianSourceDomainError("Darken factor requires 0 < loading < saturation loading")
        darken = 1.0 / (sips.pressure_exponent * (1.0 - theta))
        mtc = (
            source.diffusion_frequency_m2_s
            * np.exp(-source.activation_energy_j_mol / (R * temperature))
            * darken / PARTICLE_RADIUS_M**2
        )
        rate[index] = mtc * (target[index] ** 2 - current[index] ** 2) / (2.0 * current[index])
    return rate


def linear_valve_flow_mol_s(delta_pressure_bar: float, coefficient_kmol_s_bar: float) -> float:
    """Convert article Eq. 14 and Table 4 coefficients to engine units."""
    if coefficient_kmol_s_bar < 0.0:
        raise ValueError("Linear valve coefficient must be non-negative")
    return float(max(float(delta_pressure_bar), 0.0) * coefficient_kmol_s_bar * 1000.0)


def molecular_binary_diffusivity_o2_n2_m2_s(pressure_bar: float, temperature_k: float) -> float:
    """Chapman--Enskog O2/N2 binary diffusivity used by article Eq. 3.

    The conventional 0.0018583 expression returns cm2/s for T [K], molecular
    weight [g/mol], collision diameter [angstrom] and pressure [atm].  The
    result is converted explicitly to m2/s.  This is a dilute-gas transport
    correlation; it is not an adsorbent micropore diffusivity.
    """
    pressure = float(pressure_bar)
    temperature = float(temperature_k)
    if not np.isfinite(pressure) or pressure <= 0.0:
        raise GuderianSourceDomainError("Chapman-Enskog pressure must be finite and positive")
    _temperature_guard(temperature)
    sigma = 0.5 * (_LJ_SIGMA_ANGSTROM["O2"] + _LJ_SIGMA_ANGSTROM["N2"])
    epsilon = np.sqrt(_LJ_EPSILON_OVER_K_K["O2"] * _LJ_EPSILON_OVER_K_K["N2"])
    reduced_temperature = temperature / epsilon
    collision = (
        1.06036 / reduced_temperature**0.15610
        + 0.19300 / np.exp(0.47635 * reduced_temperature)
        + 1.03587 / np.exp(1.52996 * reduced_temperature)
        + 1.76474 / np.exp(3.89411 * reduced_temperature)
    )
    pressure_atm = pressure / 1.01325
    mass_factor = np.sqrt(
        1.0 / _MOLAR_MASS_G_MOL["O2"] + 1.0 / _MOLAR_MASS_G_MOL["N2"]
    )
    diffusivity_cm2_s = (
        0.0018583 * temperature**1.5 * mass_factor
        / (pressure_atm * sigma**2 * collision)
    )
    diffusivity = diffusivity_cm2_s * 1.0e-4
    if not np.isfinite(diffusivity) or diffusivity <= 0.0:
        raise ArithmeticError("Chapman-Enskog diffusivity is non-positive")
    return float(diffusivity)


def axial_dispersion_o2_n2_m2_s(
    superficial_velocity_m_s: float,
    pressure_bar: float,
    temperature_k: float,
    interparticle_voidage: float = 0.404,
    particle_radius_m: float = PARTICLE_RADIUS_M,
) -> np.ndarray:
    """Evaluate Marcinek et al. Eq. 3 in the fixed O2/N2 source order."""
    velocity = abs(float(superficial_velocity_m_s))
    voidage = float(interparticle_voidage)
    radius = float(particle_radius_m)
    if not np.isfinite(velocity) or velocity < 0.0:
        raise ValueError("Superficial velocity must be finite")
    if not 0.0 < voidage < 1.0 or not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("Eq. 3 requires physical voidage and particle radius")
    molecular = molecular_binary_diffusivity_o2_n2_m2_s(pressure_bar, temperature_k)
    if velocity == 0.0:
        dispersion = 0.73 * molecular
    else:
        dispersion = 0.73 * molecular + (
            velocity * radius
            / (voidage * (1.0 + 9.49 * voidage * molecular / (2.0 * velocity * radius)))
        )
    return np.full(2, float(dispersion))

