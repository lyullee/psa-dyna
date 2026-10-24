import numpy as np
import pytest

from psa_dynamic.guderian_2021_n2 import (
    adsorption_energy_magnitude_j_kg,
    adsorption_energy_temperature_derivative_j_kg_k,
    axial_dispersion_o2_n2_m2_s,
    GuderianSourceDomainError,
    heat_of_adsorption_magnitude_j_mol,
    iast_loading_mol_kg,
    linear_valve_flow_mol_s,
    molecular_binary_diffusivity_o2_n2_m2_s,
    pressure_at_spreading_pressure,
    pure_loading_mol_kg,
    qdf_rate_mol_kg_s,
    spreading_pressure,
    temperature_sips,
)


@pytest.mark.parametrize(
    "component,pressure_bar,expected_mol_kg",
    [("O2", 0.964, 0.3564), ("N2", 0.95, 0.3779)],
)
def test_reference_temperature_sips_matches_table_scale(component, pressure_bar, expected_mol_kg):
    calculated = pure_loading_mol_kg(component, pressure_bar, 293.15)
    assert calculated == pytest.approx(expected_mol_kg, rel=0.025)


def test_spreading_pressure_inverse_roundtrip():
    for name in ("O2", "N2"):
        for pressure in (0.01, 1.0, 8.0, 20.0):
            pi = spreading_pressure(name, pressure, 313.15)
            assert pressure_at_spreading_pressure(name, pi, 313.15) == pytest.approx(pressure, rel=2e-12)


def test_iast_pure_limit_and_binary_inventory():
    pure = iast_loading_mol_kg([2.0, 0.0], 293.15)
    assert pure[0] == pytest.approx(pure_loading_mol_kg("O2", 2.0, 293.15), rel=1e-12)
    assert pure[1] == 0.0
    mixture = iast_loading_mol_kg([1.6, 6.4], 293.15)
    assert mixture.shape == (2,)
    assert np.all(mixture > 0.0)
    assert mixture.sum() < sum(temperature_sips(name, 293.15).saturation_mol_kg for name in ("O2", "N2"))


def test_qdf_is_zero_at_equilibrium_and_oxygen_is_kinetically_faster():
    equilibrium = iast_loading_mol_kg([1.6, 6.4], 293.15)
    assert np.allclose(qdf_rate_mol_kg_s(equilibrium, equilibrium, 293.15), 0.0, atol=1e-14)
    current = 0.5 * equilibrium
    rate = qdf_rate_mol_kg_s(current, equilibrium, 293.15)
    assert np.all(rate > 0.0)
    assert rate[0] / equilibrium[0] > rate[1] / equilibrium[1]


def test_heat_and_linear_valve_units_are_explicit():
    assert heat_of_adsorption_magnitude_j_mol("O2", 1.0, 293.15) == pytest.approx(13067.46, rel=2e-4)
    assert heat_of_adsorption_magnitude_j_mol("N2", 1.0, 293.15) == pytest.approx(12273.35, rel=2e-4)
    assert linear_valve_flow_mol_s(0.25, 2.5e-5) == pytest.approx(0.00625)


def test_integrated_adsorption_energy_derivatives_match_heat_and_temperature_difference():
    loading = 1.0
    temperature = 303.15
    dq = 1e-5
    dt = 1e-3
    for name in ("O2", "N2"):
        numerical_q = (
            adsorption_energy_magnitude_j_kg(name, loading + dq, temperature)
            - adsorption_energy_magnitude_j_kg(name, loading - dq, temperature)
        ) / (2.0 * dq)
        assert numerical_q == pytest.approx(
            heat_of_adsorption_magnitude_j_mol(name, loading, temperature), rel=2e-8
        )
        numerical_t = (
            adsorption_energy_magnitude_j_kg(name, loading, temperature + dt)
            - adsorption_energy_magnitude_j_kg(name, loading, temperature - dt)
        ) / (2.0 * dt)
        assert numerical_t == pytest.approx(
            adsorption_energy_temperature_derivative_j_kg_k(name, loading, temperature), rel=2e-7
        )


def test_source_domain_is_fail_closed():
    with pytest.raises(GuderianSourceDomainError):
        pure_loading_mol_kg("O2", 20.1, 293.15)
    with pytest.raises(GuderianSourceDomainError):
        iast_loading_mol_kg([5.0, 16.0], 293.15)
    with pytest.raises(GuderianSourceDomainError):
        qdf_rate_mol_kg_s([0.0, 0.1], [0.2, 0.2], 293.15)


def test_chapman_enskog_binary_diffusivity_has_expected_reference_scale():
    diffusivity = molecular_binary_diffusivity_o2_n2_m2_s(1.01325, 293.15)
    assert diffusivity == pytest.approx(1.982209e-5, rel=5e-6)
    assert molecular_binary_diffusivity_o2_n2_m2_s(2.0265, 293.15) == pytest.approx(
        0.5 * diffusivity, rel=2e-12
    )


def test_source_eq3_dispersion_is_dynamic_and_binary_shared():
    molecular = molecular_binary_diffusivity_o2_n2_m2_s(1.01325, 293.15)
    stagnant = axial_dispersion_o2_n2_m2_s(0.0, 1.01325, 293.15)
    flowing = axial_dispersion_o2_n2_m2_s(0.08, 1.01325, 293.15)
    high_pressure = axial_dispersion_o2_n2_m2_s(0.08, 8.0, 293.15)
    assert stagnant.shape == (2,)
    assert np.all(stagnant == pytest.approx(0.73 * molecular, rel=2e-12))
    assert np.all(flowing > stagnant)
    # At this finite velocity Eq. (3)'s mechanical-dispersion denominator
    # decreases faster than its molecular term when pressure increases.  Do
    # not impose a guessed monotonic direction; verify the published equation
    # remains pressure-sensitive and uses one binary coefficient for O2/N2.
    assert np.all(high_pressure > flowing)
    assert flowing[0] == pytest.approx(flowing[1], rel=0.0, abs=0.0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"pressure_bar": 0.0},
        {"pressure_bar": -1.0},
        {"particle_radius_m": 0.0},
        {"interparticle_voidage": 1.0},
    ],
)
def test_source_eq3_dispersion_rejects_nonphysical_inputs(kwargs):
    arguments = {
        "superficial_velocity_m_s": 0.08,
        "pressure_bar": 1.01325,
        "temperature_k": 293.15,
    }
    arguments.update(kwargs)
    with pytest.raises((ValueError, GuderianSourceDomainError)):
        axial_dispersion_o2_n2_m2_s(**arguments)

