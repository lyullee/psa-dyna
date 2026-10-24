"""Source-neutral Peng-Robinson gas-phase primitives for PSA adapters.

All inputs are SI: pressure in Pa, temperature in K, critical pressure in Pa,
and molar density in mol/m3.  No adsorbent or process parameters live here.
The module is intentionally separate from ``core.py`` until each source-local
adapter supplies and tests its own critical-property and binary-interaction
tables.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


R = 8.31446261815324


@dataclass(frozen=True)
class PengRobinsonComponents:
    names: tuple[str, ...]
    critical_temperature_k: tuple[float, ...]
    critical_pressure_pa: tuple[float, ...]
    acentric_factor: tuple[float, ...]
    binary_interaction: tuple[tuple[float, ...], ...] | None = None

    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        tc = np.asarray(self.critical_temperature_k, dtype=float)
        pc = np.asarray(self.critical_pressure_pa, dtype=float)
        omega = np.asarray(self.acentric_factor, dtype=float)
        if len(self.names) == 0 or tc.shape != (len(self.names),) or pc.shape != tc.shape or omega.shape != tc.shape:
            raise ValueError("Peng-Robinson component vectors must have one common nonzero length")
        if np.any(~np.isfinite(tc)) or np.any(~np.isfinite(pc)) or np.any(~np.isfinite(omega)):
            raise ValueError("Peng-Robinson critical properties must be finite")
        if np.any(tc <= 0.0) or np.any(pc <= 0.0) or len(set(self.names)) != len(self.names):
            raise ValueError("Peng-Robinson requires positive critical properties and unique names")
        if self.binary_interaction is None:
            kij = np.zeros((len(self.names), len(self.names)), dtype=float)
        else:
            kij = np.asarray(self.binary_interaction, dtype=float)
            if kij.shape != (len(self.names), len(self.names)) or np.any(~np.isfinite(kij)):
                raise ValueError("Peng-Robinson binary-interaction matrix has the wrong shape")
            if not np.allclose(kij, kij.T, rtol=0.0, atol=1e-14) or not np.allclose(np.diag(kij), 0.0, rtol=0.0, atol=1e-14):
                raise ValueError("Peng-Robinson binary-interaction matrix must be symmetric with a zero diagonal")
        return tc, pc, omega, kij


def _composition(composition: np.ndarray, count: int) -> np.ndarray:
    y = np.asarray(composition)
    if y.shape != (count,) or np.any(~np.isfinite(np.real(y))) or np.any(np.real(y) < 0.0):
        raise ValueError("Peng-Robinson composition must be a finite nonnegative component vector")
    total = np.sum(y)
    if np.real(total) <= 0.0:
        raise ValueError("Peng-Robinson composition sum must be positive")
    return y / total


def mixture_parameters(
    temperature_k: float | complex,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> tuple[complex | float, complex | float]:
    """Return ``a_alpha`` [Pa m6/mol2] and ``b`` [m3/mol]."""
    tc, pc, omega, kij = components.arrays()
    if np.real(temperature_k) <= 0.0 or not np.isfinite(np.real(temperature_k)):
        raise ValueError("Peng-Robinson temperature must be finite and positive")
    y = _composition(composition, len(tc))
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1.0 + kappa * (1.0 - np.sqrt(temperature_k / tc))) ** 2
    ai = 0.45724 * R**2 * tc**2 / pc * alpha
    bi = 0.07780 * R * tc / pc
    cross = np.sqrt(ai[:, None] * ai[None, :]) * (1.0 - kij)
    return np.sum(y[:, None] * y[None, :] * cross), np.dot(y, bi)


def mixture_parameters_with_temperature_derivative(
    temperature_k: float | complex,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> tuple[complex | float, complex | float, complex | float]:
    """Return PR ``a_alpha``, ``b`` and ``d(a_alpha)/dT`` at fixed composition."""
    tc, pc, omega, kij = components.arrays()
    if np.real(temperature_k) <= 0.0 or not np.isfinite(np.real(temperature_k)):
        raise ValueError("Peng-Robinson temperature must be finite and positive")
    y = _composition(composition, len(tc))
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    factor = 1.0 + kappa * (1.0 - np.sqrt(temperature_k / tc))
    a0 = 0.45724 * R**2 * tc**2 / pc
    ai = a0 * factor**2
    dai_dt = -a0 * kappa * factor / np.sqrt(temperature_k * tc)
    bi = 0.07780 * R * tc / pc
    cross = np.sqrt(ai[:, None] * ai[None, :]) * (1.0 - kij)
    cross_derivative = 0.5 * cross * (
        dai_dt[:, None] / ai[:, None] + dai_dt[None, :] / ai[None, :]
    )
    weights = y[:, None] * y[None, :]
    return np.sum(weights * cross), np.dot(y, bi), np.sum(weights * cross_derivative)


def molar_internal_energy_departure_from_density(
    molar_density_mol_m3: float | complex,
    temperature_k: float | complex,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> complex | float:
    """Return the PR residual molar internal energy [J/mol].

    The expression is the standard cubic-EOS departure integral at fixed
    composition.  It tends continuously to zero as density tends to zero.
    """
    if np.real(molar_density_mol_m3) <= 0.0 or not np.isfinite(np.real(molar_density_mol_m3)):
        raise ValueError("Peng-Robinson molar density must be finite and positive")
    a_mix, b_mix, da_dt = mixture_parameters_with_temperature_derivative(
        temperature_k, composition, components
    )
    molar_volume = 1.0 / molar_density_mol_m3
    if np.real(molar_volume - b_mix) <= 0.0:
        raise ValueError("Peng-Robinson molar volume must exceed the co-volume")
    root_two = np.sqrt(2.0)
    # log1p form avoids subtractive loss in the dilute-gas limit.
    logarithm = np.log1p(
        2.0 * root_two * b_mix
        / (molar_volume + (1.0 - root_two) * b_mix)
    )
    return (temperature_k * da_dt - a_mix) / (2.0 * root_two * b_mix) * logarithm


def molar_enthalpy_departure(
    pressure_pa: float,
    temperature_k: float,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> float:
    """Return PR residual molar enthalpy [J/mol] on the gas root."""
    density = molar_density_from_pressure(pressure_pa, temperature_k, composition, components)
    z_factor = compressibility_factor(pressure_pa, temperature_k, composition, components)
    internal = molar_internal_energy_departure_from_density(
        density, temperature_k, composition, components
    )
    return float(np.real(internal + (z_factor - 1.0) * R * temperature_k))


def total_internal_energy_departure(
    inventory_mol: np.ndarray,
    volume_m3: float,
    temperature_k: float | complex,
    components: PengRobinsonComponents,
) -> complex | float:
    """Return total PR residual internal energy [J] for a rigid gas volume."""
    inventory = np.asarray(inventory_mol)
    if inventory.shape != (len(components.names),):
        raise ValueError("Peng-Robinson inventory vector has the wrong shape")
    total = np.sum(inventory)
    if np.real(total) <= 0.0 or not np.isfinite(np.real(total)):
        raise ValueError("Peng-Robinson total inventory must be finite and positive")
    if not np.isfinite(volume_m3) or volume_m3 <= 0.0:
        raise ValueError("Peng-Robinson gas volume must be finite and positive")
    composition = inventory / total
    return total * molar_internal_energy_departure_from_density(
        total / volume_m3, temperature_k, composition, components
    )


def total_internal_energy_departure_derivatives(
    inventory_mol: np.ndarray,
    volume_m3: float,
    temperature_k: float,
    components: PengRobinsonComponents,
) -> tuple[np.ndarray, float]:
    """Analytic derivatives of total residual energy with respect to ``n,T``.

    At fixed rigid volume the conventional PR residual internal energy can be
    rearranged into a particularly compact inventory form::

        U^R = S/(2 sqrt(2) B) ln[(V+(1+sqrt(2))B)/(V+(1-sqrt(2))B)]

    where ``B = sum(n_i b_i)`` and
    ``S = n.T @ (T*dA/dT - A) @ n``.  Differentiating this expression avoids
    one complex-step EOS evaluation per species and per spatial cell at every
    ODE right-hand-side call.  The returned quantities remain the exact
    partial derivatives required by the conservative internal-energy balance.
    """
    inventory = np.asarray(inventory_mol, dtype=float)
    if inventory.shape != (len(components.names),) or np.any(inventory <= 0.0):
        raise ValueError("PR departure derivatives require positive real inventories")
    if not np.isfinite(volume_m3) or volume_m3 <= 0.0:
        raise ValueError("Peng-Robinson gas volume must be finite and positive")
    if not np.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("Peng-Robinson temperature must be finite and positive")

    tc, pc, omega, kij = components.arrays()
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    root_temperature_ratio = np.sqrt(temperature_k / tc)
    factor = 1.0 + kappa * (1.0 - root_temperature_ratio)
    a0 = 0.45724 * R**2 * tc**2 / pc
    ai = a0 * factor**2
    dai_dt = -a0 * kappa * factor / np.sqrt(temperature_k * tc)
    d2ai_dt2 = 0.5 * a0 * kappa * (
        kappa / (temperature_k * tc)
        + factor / (temperature_k * np.sqrt(temperature_k * tc))
    )
    bi = 0.07780 * R * tc / pc

    attraction = np.sqrt(ai[:, None] * ai[None, :]) * (1.0 - kij)
    logarithmic_slope = dai_dt / ai
    attraction_dt = 0.5 * attraction * (
        logarithmic_slope[:, None] + logarithmic_slope[None, :]
    )
    logarithmic_curvature = d2ai_dt2 / ai - logarithmic_slope**2
    attraction_dt2 = attraction * (
        0.25 * (logarithmic_slope[:, None] + logarithmic_slope[None, :]) ** 2
        + 0.5 * (
            logarithmic_curvature[:, None] + logarithmic_curvature[None, :]
        )
    )

    b_inventory = float(np.dot(inventory, bi))
    if b_inventory <= 0.0:
        raise FloatingPointError("PR inventory co-volume must be positive")
    energy_matrix = temperature_k * attraction_dt - attraction
    energy_quadratic = float(inventory @ energy_matrix @ inventory)
    root_two = np.sqrt(2.0)
    numerator = volume_m3 + (1.0 + root_two) * b_inventory
    denominator = volume_m3 + (1.0 - root_two) * b_inventory
    if numerator <= 0.0 or denominator <= 0.0:
        raise FloatingPointError("PR departure-energy logarithm is outside its domain")
    logarithm = np.log(numerator / denominator)
    logarithm_derivative = (
        (1.0 + root_two) / numerator
        - (1.0 - root_two) / denominator
    )
    prefactor = 1.0 / (2.0 * root_two)
    quadratic_gradient = 2.0 * (energy_matrix @ inventory)
    dn = prefactor * (
        logarithm
        * (quadratic_gradient / b_inventory
           - energy_quadratic * bi / b_inventory**2)
        + energy_quadratic / b_inventory * logarithm_derivative * bi
    )
    energy_matrix_dt = temperature_k * attraction_dt2
    dt = prefactor * logarithm / b_inventory * float(
        inventory @ energy_matrix_dt @ inventory
    )
    if np.any(~np.isfinite(dn)) or not np.isfinite(dt):
        raise FloatingPointError("PR departure-energy derivatives are non-finite")
    return dn, float(dt)


def pressure_from_molar_density(
    molar_density_mol_m3: float | complex,
    temperature_k: float | complex,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> complex | float:
    """Evaluate pressure directly from the conservative molar inventory."""
    if np.real(molar_density_mol_m3) <= 0.0 or not np.isfinite(np.real(molar_density_mol_m3)):
        raise ValueError("Peng-Robinson molar density must be finite and positive")
    a_mix, b_mix = mixture_parameters(temperature_k, composition, components)
    molar_volume = 1.0 / molar_density_mol_m3
    if np.real(molar_volume - b_mix) <= 0.0:
        raise ValueError("Peng-Robinson molar volume must exceed the co-volume")
    denominator = molar_volume * (molar_volume + b_mix) + b_mix * (molar_volume - b_mix)
    return R * temperature_k / (molar_volume - b_mix) - a_mix / denominator


def compressibility_factor(
    pressure_pa: float,
    temperature_k: float,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> float:
    """Return the largest admissible real gas root of the PR cubic."""
    if not np.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise ValueError("Peng-Robinson pressure must be finite and positive")
    a_mix, b_mix = mixture_parameters(temperature_k, composition, components)
    a_mix = float(np.real(a_mix))
    b_mix = float(np.real(b_mix))
    A = a_mix * pressure_pa / (R * temperature_k) ** 2
    B = b_mix * pressure_pa / (R * temperature_k)
    roots = np.roots([1.0, -(1.0 - B), A - 3.0 * B**2 - 2.0 * B, -(A * B - B**2 - B**3)])
    candidates = sorted(float(root.real) for root in roots if abs(root.imag) <= 1e-10 and root.real > B)
    if not candidates:
        raise FloatingPointError("Peng-Robinson cubic has no admissible real gas root")
    return candidates[-1]


def molar_density_from_pressure(
    pressure_pa: float,
    temperature_k: float,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> float:
    z_factor = compressibility_factor(pressure_pa, temperature_k, composition, components)
    return pressure_pa / (z_factor * R * temperature_k)


def cubic_residual(
    z_factor: float,
    pressure_pa: float,
    temperature_k: float,
    composition: np.ndarray,
    components: PengRobinsonComponents,
) -> float:
    a_mix, b_mix = mixture_parameters(temperature_k, composition, components)
    A = float(np.real(a_mix)) * pressure_pa / (R * temperature_k) ** 2
    B = float(np.real(b_mix)) * pressure_pa / (R * temperature_k)
    return float(
        z_factor**3
        - (1.0 - B) * z_factor**2
        + (A - 3.0 * B**2 - 2.0 * B) * z_factor
        - (A * B - B**2 - B**3)
    )
