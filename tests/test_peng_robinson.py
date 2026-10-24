import numpy as np
import pytest

from psa_dynamic.core import Bed
from psa_dynamic.peng_robinson import (
    PengRobinsonComponents,
    compressibility_factor,
    cubic_residual,
    molar_enthalpy_departure,
    molar_density_from_pressure,
    molar_internal_energy_departure_from_density,
    pressure_from_molar_density,
    total_internal_energy_departure,
    total_internal_energy_departure_derivatives,
)


AIR = PengRobinsonComponents(
    names=("O2", "N2"),
    critical_temperature_k=(154.581, 126.192),
    critical_pressure_pa=(5.043e6, 3.3958e6),
    acentric_factor=(0.0222, 0.0372),
)


@pytest.mark.parametrize("pressure_bar", [0.1, 1.0, 6.0, 8.0, 10.0, 20.0])
def test_pr_air_cubic_and_density_roundtrip(pressure_bar):
    pressure = pressure_bar * 1e5
    composition = np.array([0.209, 0.791])
    z_factor = compressibility_factor(pressure, 293.15, composition, AIR)
    assert abs(cubic_residual(z_factor, pressure, 293.15, composition, AIR)) < 1e-11
    density = molar_density_from_pressure(pressure, 293.15, composition, AIR)
    recovered = pressure_from_molar_density(density, 293.15, composition, AIR)
    assert recovered == pytest.approx(pressure, rel=2e-12)
    assert 0.97 < z_factor < 1.01


def test_pr_low_pressure_approaches_ideal_gas():
    composition = np.array([0.209, 0.791])
    z_factor = compressibility_factor(1000.0, 293.15, composition, AIR)
    assert z_factor == pytest.approx(1.0, abs=5e-5)


def test_pr_departure_energy_obeys_enthalpy_internal_energy_identity():
    pressure = 8.0e5
    temperature = 293.15
    composition = np.array([0.209, 0.791])
    density = molar_density_from_pressure(pressure, temperature, composition, AIR)
    internal = molar_internal_energy_departure_from_density(
        density, temperature, composition, AIR
    )
    enthalpy = molar_enthalpy_departure(pressure, temperature, composition, AIR)
    z_factor = compressibility_factor(pressure, temperature, composition, AIR)
    assert enthalpy - internal == pytest.approx(
        (z_factor - 1.0) * 8.31446261815324 * temperature, rel=2e-12, abs=2e-10
    )
    assert internal < 0.0


def test_pr_total_departure_derivatives_match_centered_finite_differences():
    inventory = np.array([0.15, 0.55])
    volume = 0.002
    temperature = 303.15
    dn, dt = total_internal_energy_departure_derivatives(
        inventory, volume, temperature, AIR
    )
    for index in range(2):
        step = 1e-6 * inventory[index]
        plus = inventory.copy(); plus[index] += step
        minus = inventory.copy(); minus[index] -= step
        numerical = (
            total_internal_energy_departure(plus, volume, temperature, AIR)
            - total_internal_energy_departure(minus, volume, temperature, AIR)
        ) / (2.0 * step)
        assert dn[index] == pytest.approx(numerical, rel=2e-7, abs=2e-6)
    step_t = 1e-4
    numerical_t = (
        total_internal_energy_departure(inventory, volume, temperature + step_t, AIR)
        - total_internal_energy_departure(inventory, volume, temperature - step_t, AIR)
    ) / (2.0 * step_t)
    assert dt == pytest.approx(numerical_t, rel=2e-7, abs=2e-6)


@pytest.mark.parametrize(
    "inventory,volume,temperature",
    [
        ([0.004, 0.019], 8.0e-4, 278.15),
        ([0.31, 0.07], 3.0e-3, 293.15),
        ([0.12, 0.88], 8.0e-3, 333.15),
    ],
)
def test_pr_analytic_departure_derivatives_cover_source_temperature_and_composition_range(
    inventory, volume, temperature
):
    """Guard the analytic speed-up over the complete documented fit range."""
    inventory = np.asarray(inventory, dtype=float)
    dn, dt = total_internal_energy_departure_derivatives(
        inventory, volume, temperature, AIR
    )
    for index in range(2):
        step = max(1.0e-9, 2.0e-6 * inventory[index])
        plus = inventory.copy(); plus[index] += step
        minus = inventory.copy(); minus[index] -= step
        numerical = (
            total_internal_energy_departure(plus, volume, temperature, AIR)
            - total_internal_energy_departure(minus, volume, temperature, AIR)
        ) / (2.0 * step)
        assert dn[index] == pytest.approx(numerical, rel=8e-7, abs=1e-5)
    step_t = 2.0e-4
    numerical_t = (
        total_internal_energy_departure(inventory, volume, temperature + step_t, AIR)
        - total_internal_energy_departure(inventory, volume, temperature - step_t, AIR)
    ) / (2.0 * step_t)
    assert dt == pytest.approx(numerical_t, rel=8e-7, abs=1e-5)


def test_pr_rejects_unqualified_component_or_mixture_data():
    with pytest.raises(ValueError):
        compressibility_factor(8e5, 293.15, np.array([1.0]), AIR)
    invalid = PengRobinsonComponents(
        names=("O2", "N2"),
        critical_temperature_k=(154.581, 126.192),
        critical_pressure_pa=(5.043e6, 3.3958e6),
        acentric_factor=(0.0222, 0.0372),
        binary_interaction=((0.0, 0.1), (0.2, 0.0)),
    )
    with pytest.raises(ValueError):
        compressibility_factor(8e5, 293.15, np.array([0.209, 0.791]), invalid)


def _bed(equation_of_state):
    species = [
        {"name": "O2", "mw_kg_mol": 0.031998, "cp_j_mol_k": 29.38,
         "critical_temperature_k": 154.581, "critical_pressure_pa": 5.043e6, "acentric_factor": 0.0222},
        {"name": "N2", "mw_kg_mol": 0.0280134, "cp_j_mol_k": 29.12,
         "critical_temperature_k": 126.192, "critical_pressure_pa": 3.3958e6, "acentric_factor": 0.0372},
    ]
    material = {
        "end_fraction": 1.0,
        "a1": [[1.0], [1.0]],
        "a2": [[0.0], [0.0]],
        "b0_per_bar": [[0.1], [0.1]],
        "affinity_energy_j_mol": [[0.0], [0.0]],
        "binding_energy_j_mol": [[10000.0], [10000.0]],
        "ldf_per_s": [0.1, 0.1],
    }
    spec = {
        "n_cells": 2,
        "length_m": 0.58,
        "diameter_m": 0.066,
        "void_fraction": 0.404,
        "particle_diameter_m": 0.00166,
        "bulk_density_kg_m3": 711.0,
        "solid_cp_j_kg_k": 880.0,
        "ambient_k": 293.15,
        "materials": [material],
        "equation_of_state": equation_of_state,
    }
    return Bed(spec, species)


def test_bed_pr_inventory_roundtrips_reference_pressure_and_preserves_ideal_default():
    composition = np.array([0.209, 0.791])
    pr_bed = _bed("peng_robinson")
    ideal_bed = _bed("ideal_gas")
    pr_state = pr_bed.initial(8.0, composition, 293.15)
    ideal_state = ideal_bed.initial(8.0, composition, 293.15)
    assert np.allclose(pr_bed.fields(pr_state)[4], 8.0, rtol=0.0, atol=2e-12)
    assert np.allclose(ideal_bed.fields(ideal_state)[4], 8.0, rtol=0.0, atol=2e-12)
    assert pr_state.reshape(2, -1)[:, :2].sum() > ideal_state.reshape(2, -1)[:, :2].sum()


def test_bed_pr_scalar_and_vectorized_ergun_flux_are_equivalent():
    composition = np.array([0.209, 0.791])
    bed = _bed("peng_robinson")
    state = bed.initial(8.0, composition, 293.15).reshape(2, -1)
    density_7bar = molar_density_from_pressure(7e5, 293.15, composition, AIR)
    state[1, :2] = composition * density_7bar * bed.vgas[1]
    scalar_flux, scalar_heat = bed.internal_flux_scalar(state.ravel())
    vector_flux, vector_heat = bed.internal_flux(state.ravel())
    assert np.allclose(vector_flux, scalar_flux, rtol=2e-13, atol=1e-13)
    assert np.allclose(vector_heat, scalar_heat, rtol=2e-13, atol=1e-10)
