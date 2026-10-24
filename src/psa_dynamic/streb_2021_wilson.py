"""Explicit-matrix, *unqualified* Wilson/RAST equilibrium calculation.

The caller must supply every Wilson Lambda_ij explicitly.  Streb and
Mazzotti (2021, part 2) report a fitted CO2/CH4 interaction energy, but the
source contract does not yet establish the nu_j/nu_i ratio in Eq. (2).
Consequently this module does not infer a source matrix. Its result can be
used by the loading-LDF engine as a diagnostic candidate, but never as a
source-validated design or training label.

Equilibrium follows y_i P = x_i gamma_i P_i^0 at a common reduced spreading
integral.  For a Wilson excess Gibbs energy independent of spreading pressure,
the excess contribution to 1/q_total is zero; the harmonic pure-loading sum
is then valid.  Pressure dependence of gamma would need an additional term.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.special import softmax

from .streb_2021_sips import (
    inverse_pure_spreading_pressure_bar,
    pure_loading_mol_kg,
    pure_spreading_integral_mol_kg,
)


SOURCE_DOI = "https://doi.org/10.1007/s10450-021-00308-w"
STATUS = "mathematical_candidate_not_source_validated"


def _validated_matrix(matrix: Sequence[Sequence[float]], count: int) -> np.ndarray:
    result = np.asarray(matrix, dtype=float)
    if result.shape != (count, count):
        raise ValueError("Wilson Lambda matrix shape must match the species order")
    if not np.all(np.isfinite(result)) or np.any(result <= 0):
        raise ValueError("All Wilson Lambda entries must be finite and positive")
    if not np.allclose(np.diag(result), 1.0, rtol=0, atol=1e-12):
        raise ValueError("Wilson Lambda diagonal must be one")
    return result


def wilson_log_activity_coefficients(
    adsorbed_fractions: Sequence[float], matrix: Sequence[Sequence[float]]
) -> np.ndarray:
    """Evaluate ln(gamma_i) using the conventional Wilson row convention.

    Lambda[i,j] multiplies x[j] in row sum i.  A matrix of ones is exactly
    the IAST limit; the pure-component activity coefficient equals one.
    """
    x = np.asarray(adsorbed_fractions, dtype=float)
    if x.ndim != 1 or len(x) == 0 or not np.all(np.isfinite(x)):
        raise ValueError("Adsorbed fractions must be a finite one-dimensional vector")
    if np.any(x < 0) or not math.isclose(float(np.sum(x)), 1.0, rel_tol=0, abs_tol=1e-10):
        raise ValueError("Adsorbed fractions must be nonnegative and sum to one")
    lam = _validated_matrix(matrix, len(x))
    return _wilson_log_activity_raw(x, lam)


def _wilson_log_activity_raw(x: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """Validated-input inner kernel, avoiding matrix checks in each root call."""
    row_sums = lam @ x
    # The second term for species i is sum_j x_j Lambda[j,i] / row_sum[j].
    return 1.0 - np.log(row_sums) - lam.T @ (x / row_sums)


def _wilson_log_activity_dx(x: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """Return d ln(gamma_i)/d x_m for the conventional Wilson rows.

    The derivative is evaluated before the softmax chain rule.  It is not
    an assumed equilibrium d q/dp or d q/dT for the bed equations.
    """
    row_sums = lam @ x
    return (-lam / row_sums[:, None]
            - lam.T / row_sums[None, :]
            + lam.T @ ((x / row_sums**2)[:, None] * lam))


def solve_wilson_rast_mixture(
    species: Sequence[str],
    gas_fractions: Sequence[float],
    pressure_bar: float,
    temperature_k: float,
    matrix: Sequence[Sequence[float]],
    *,
    parameter_provenance: str,
) -> dict:
    """Solve a source-inspired, explicit-parameter equilibrium candidate.

    No fitted or guessed interaction parameter is supplied by this function.
    Positive feed components are solved; zero-feed components receive zero
    loading and do not enter the nonlinear solve.  A poor numerical root
    fails closed instead of silently returning an approximate equilibrium.
    """
    names = tuple(species)
    y = np.asarray(gas_fractions, dtype=float)
    if len(names) == 0 or len(names) != len(set(names)) or set(names) - {"CO2", "CH4", "H2"}:
        raise ValueError("Provide unique species from the source Table 2")
    if y.shape != (len(names),) or not np.all(np.isfinite(y)) or np.any(y < 0):
        raise ValueError("Gas fractions must be finite and nonnegative")
    if not math.isclose(float(np.sum(y)), 1.0, rel_tol=0, abs_tol=1e-10):
        raise ValueError("Gas fractions must sum to one")
    if not math.isfinite(pressure_bar) or pressure_bar <= 0:
        raise ValueError("Pressure must be finite and positive in bar")
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("Temperature must be finite and positive in kelvin")
    if not isinstance(parameter_provenance, str) or not parameter_provenance.strip():
        raise ValueError("Explicit Wilson-matrix provenance is required")
    lam_all = _validated_matrix(matrix, len(names))
    active = np.flatnonzero(y > 0)
    active_names = tuple(names[i] for i in active)
    active_y = y[active]
    lam = lam_all[np.ix_(active, active)]
    if len(active) == 1:
        loading = pure_loading_mol_kg(active_names[0], pressure_bar, temperature_k)
        psi = pure_spreading_integral_mol_kg(active_names[0], pressure_bar, temperature_k)
        x = np.zeros(len(names))
        q = np.zeros(len(names))
        x[active[0]], q[active[0]] = 1.0, loading
        return _result(names, y, x, q, psi, pressure_bar, temperature_k,
                       lam_all, parameter_provenance, 0.0)

    def p0(name: str, psi: float) -> float:
        return inverse_pure_spreading_pressure_bar(name, psi, temperature_k)

    def ideal_sum(psi: float) -> float:
        return sum(float(active_y[i]) * pressure_bar / p0(name, psi)
                   for i, name in enumerate(active_names)) - 1.0

    lo = min(pure_spreading_integral_mol_kg(name, float(active_y[i]) * pressure_bar,
                                            temperature_k)
             for i, name in enumerate(active_names))
    hi = max(pure_spreading_integral_mol_kg(name, pressure_bar, temperature_k)
             for name in active_names)
    psi_ideal = float(brentq(ideal_sum, lo, hi, xtol=1e-13, rtol=1e-13))
    x_ideal = np.array([float(active_y[i]) * pressure_bar / p0(name, psi_ideal)
                        for i, name in enumerate(active_names)])
    x_ideal /= np.sum(x_ideal)
    ref = len(active) - 1
    z0 = np.r_[math.log(psi_ideal), np.log(x_ideal[:-1] / x_ideal[ref])]

    def residual(z: np.ndarray) -> np.ndarray:
        psi = math.exp(float(z[0]))
        x_active = softmax(np.r_[z[1:], 0.0])
        ln_gamma = _wilson_log_activity_raw(x_active, lam)
        log_p0 = np.array([math.log(p0(name, psi)) for name in active_names])
        return np.log(x_active) + ln_gamma + log_p0 - np.log(active_y * pressure_bar)

    def jacobian(z: np.ndarray) -> np.ndarray:
        """Exact residual Jacobian in log-spreading and softmax coordinates."""
        psi = math.exp(float(z[0]))
        x_active = softmax(np.r_[z[1:], 0.0])
        p0_values = np.array([p0(name, psi) for name in active_names])
        pure_loadings = np.array([pure_loading_mol_kg(name, value, temperature_k)
                                  for name, value in zip(active_names, p0_values)])
        if np.any(pure_loadings <= 0):
            raise ArithmeticError("Wilson/RAST analytic Jacobian requires positive pure loadings")
        # d psi / d ln(P_i^0) = q_i^0, hence d ln(P_i^0)/d ln(psi) = psi/q_i^0.
        first = psi / pure_loadings
        # z_k = ln(x_k/x_ref), k < ref; softmax supplies the simplex tangent.
        dx_dz = x_active[:, None] * (np.eye(len(active))[:, :ref]
                                     - x_active[None, :ref])
        dlogx_dz = np.eye(len(active))[:, :ref] - x_active[None, :ref]
        composition = dlogx_dz + _wilson_log_activity_dx(x_active, lam) @ dx_dz
        return np.column_stack((first, composition))

    lower = np.r_[max(math.log(psi_ideal) - 15.0, -50.0), np.full(ref, -40.0)]
    upper = np.r_[min(math.log(psi_ideal) + 15.0, 20.0), np.full(ref, 40.0)]
    root = least_squares(residual, z0, jac=jacobian, bounds=(lower, upper),
                         xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=300)
    maximum_log_residual = float(np.max(np.abs(residual(root.x))))
    if not root.success or maximum_log_residual > 1e-9:
        raise ArithmeticError("Wilson/RAST equilibrium did not close to 1e-9 log residual")
    psi = math.exp(float(root.x[0]))
    x_active = softmax(np.r_[root.x[1:], 0.0])
    pure_loadings = np.array([pure_loading_mol_kg(name, p0(name, psi), temperature_k)
                              for name in active_names])
    total = 1.0 / float(np.sum(x_active / pure_loadings))
    x = np.zeros(len(names))
    q = np.zeros(len(names))
    x[active] = x_active
    q[active] = total * x_active
    if not np.all(np.isfinite(q)) or np.any(q < 0):
        raise ArithmeticError("Wilson/RAST equilibrium produced invalid loadings")
    return _result(names, y, x, q, psi, pressure_bar, temperature_k,
                   lam_all, parameter_provenance, maximum_log_residual)


def _result(names, y, x, q, psi, pressure_bar, temperature_k, lam,
            provenance, max_residual) -> dict:
    return {
        "status": STATUS,
        "species_order": list(names),
        "gas_mole_fractions": y.tolist(),
        "adsorbed_mole_fractions": x.tolist(),
        "component_loadings_mol_kg": q.tolist(),
        "total_loading_mol_kg": float(np.sum(q)),
        "reduced_spreading_integral_mol_kg": float(psi),
        "pressure_bar": float(pressure_bar),
        "temperature_k": float(temperature_k),
        "wilson_lambda_matrix": lam.tolist(),
        "parameter_provenance": provenance,
        "max_log_equilibrium_residual": float(max_residual),
        "source_doi": SOURCE_DOI,
        "release_flags": {
            "source_validation_released": False,
            "design_eligible": False,
            "dwsim_application_allowed": False,
            "tea_lca_allowed": False,
            "training_eligible": False,
        },
    }
