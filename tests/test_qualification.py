import numpy as np
import pytest
from psa_dynamic.cases import s015_case,wang_case
from psa_dynamic.core import equilibrium
from psa_dynamic.vectorized import equilibrium_batch
from psa_dynamic.qualification import MaterialDomainError,case_screen


def test_mathematical_capacity_limit_is_not_mislabelled_experimental_domain():
    case=s015_case(n_cells=4)
    audit=case_screen(case)
    co=[r for r in audit["materials"][1]["mathematical_positive_capacity_limits"] if r["species"]=="CO"][0]
    assert co["zero_capacity_temperature_k"]==pytest.approx(378.447284345)
    assert not audit["materials"][1]["experimental_fit_domain_documented"]
    assert not audit["engineering_evidence_complete"]
    assert not audit["training_eligible"]


def test_material_rejection_names_material_species_and_local_cell():
    m=s015_case()["bed"]["materials"][1]
    with pytest.raises(MaterialDomainError) as scalar:
        equilibrium(np.ones(5),379.,m)
    assert scalar.value.report["species_index"]==3
    assert scalar.value.report["temperature_k"]==379.
    with pytest.raises(MaterialDomainError) as batch:
        equilibrium_batch(np.ones((3,5)),np.array([310.,320.,379.]),m)
    assert batch.value.report["local_material_cell_index"]==2
    assert batch.value.report["species_index"]==3
    assert batch.value.report["capacity_mol_kg"]<0.


def test_flow_screen_is_a_necessary_not_sufficient_instantaneous_condition():
    case=s015_case(n_cells=4)
    original=case_screen(case)["hydraulic_necessary_screen"]
    assert original["optimistic_feed_limit_during_product_export_mol_s"]==pytest.approx(14.91134519)
    assert not original["instantaneous_reference_cap_passes_necessary_screen"]
    case["valves"]["feed"]["area_cd_m2"]*=2
    enlarged=case_screen(case)
    assert enlarged["hydraulic_necessary_screen"]["instantaneous_reference_cap_passes_necessary_screen"]
    # A bigger illustrative valve cannot qualify a plant or release training
    # data. No installed hardware or upstream-average matching was supplied.
    assert not enlarged["engineering_evidence_complete"]
    assert not enlarged["training_eligible"]


def test_material_domain_requires_grade_model_and_exact_species_coverage():
    case=s015_case(n_cells=4)
    names=[item["name"] for item in case["species"]]
    case["bed"]["materials"][0]["material_grade"]="synthetic_AC_test_grade"
    case["bed"]["materials"][0]["experimental_fit_domain"]={
        "temperature_k":[293.,323.], "pressure_bar":[0.,20.],
        "source":"synthetic test fixture", "locator":"test_qualification.py",
        "material_grade":"synthetic_AC_test_grade", "components":names,
        "equilibrium_model":"dual_site_langmuir"
    }
    audit=case_screen(case)
    assert audit["materials"][0]["experimental_fit_domain_documented"]
    assert audit["materials"][0]["experimental_fit_domain_component_match"]
    assert audit["materials"][0]["experimental_fit_domain_grade_match"]
    case["bed"]["materials"][0]["experimental_fit_domain"]["material_grade"]="unrelated_5A_grade"
    wrong_grade=case_screen(case)
    assert not wrong_grade["materials"][0]["experimental_fit_domain_documented"]
    assert not wrong_grade["materials"][0]["experimental_fit_domain_grade_match"]
    case["bed"]["materials"][0]["experimental_fit_domain"]["material_grade"]="synthetic_AC_test_grade"
    case["bed"]["materials"][0]["experimental_fit_domain"]["components"]=["H2"]
    rejected=case_screen(case)
    assert not rejected["materials"][0]["experimental_fit_domain_documented"]
    assert not rejected["materials"][0]["experimental_fit_domain_component_match"]
    case["bed"]["materials"][0]["experimental_fit_domain"]["components"]=names
    case["bed"]["materials"][0]["experimental_fit_domain"]["equilibrium_model"]="Sips"
    wrong_model=case_screen(case)
    assert not wrong_model["materials"][0]["experimental_fit_domain_equilibrium_model_match"]
    assert not wrong_model["materials"][0]["experimental_fit_domain_documented"]


def test_documented_material_domain_rejects_operating_pressure_extrapolation():
    case = s015_case("skarstrom_4step", 4)
    names = [item["name"] for item in case["species"]]
    case["bed"]["materials"][0]["material_grade"] = "synthetic_AC_test_grade"
    case["bed"]["materials"][0]["experimental_fit_domain"] = {
        "temperature_k": [293.0, 323.0],
        "pressure_bar": [0.0, 10.0],
        "source": "synthetic test fixture",
        "locator": "test_qualification.py",
        "material_grade": "synthetic_AC_test_grade",
        "components": names,
        "equilibrium_model": "dual_site_langmuir",
    }
    audit = case_screen(case)
    row = audit["materials"][0]
    assert row["experimental_fit_domain_documented"]
    assert not row["operating_point_within_documented_domain"]
    assert "declared_operating_point_outside_documented_material_domain" in audit["open_engineering_evidence_gates"]


def test_wang_domain_uses_source_component_partial_pressure_ranges():
    case = wang_case("equalized_8step", 2)
    audit = case_screen(case)
    row = audit["materials"][0]
    assert row["experimental_fit_domain_documented"]
    assert row["experimental_fit_domain"]["component_partial_pressure_bar"] == {
        "H2": [0.0, 7.0], "CH4": [0.0, 12.0],
        "N2": [0.0, 1.2], "CO2": [0.0, 1.2],
    }
    assert row["operating_point_within_documented_domain"]
    assert "experimental_isotherm_temperature_pressure_and_material_grade_domain" not in audit["open_engineering_evidence_gates"]


def test_partial_pressure_domain_requires_exact_component_coverage():
    case = wang_case("skarstrom_4step", 2)
    del case["bed"]["materials"][0]["experimental_fit_domain"]["component_partial_pressure_bar"]["CO2"]
    audit = case_screen(case)
    assert not audit["materials"][0]["experimental_fit_domain_documented"]
    assert "experimental_isotherm_temperature_pressure_and_material_grade_domain" in audit["open_engineering_evidence_gates"]


def test_inventory_and_bed_count_are_reported_before_cycle_comparisons():
    four=case_screen(s015_case("skarstrom_4step",4))
    eight=case_screen(s015_case("equalized_8step",4))
    assert (four["number_of_beds"],eight["number_of_beds"])==(2,4)
    assert eight["total_adsorbent_mass_kg"]==pytest.approx(2*four["total_adsorbent_mass_kg"])
    assert case_screen(wang_case())["hydraulic_necessary_screen"]["mode"]=="wang_controller"


def test_s015_geometry_heuristic_cannot_count_as_plant_thermal_validation():
    case=s015_case("skarstrom_4step",4)
    audit=case_screen(case)
    thermal=audit["thermal_geometry_screen"]
    assert thermal["case_diameter_m"]==pytest.approx(1.5)
    assert thermal["reference_lab_diameter_m"]==pytest.approx(0.02118)
    assert thermal["case_wall_coefficient_w_m2_k"]==pytest.approx(0.32)
    assert "plant_scale_radial_heat_transfer_validation" in audit["open_engineering_evidence_gates"]
    case["bed"]["h_bed_wall_w_m2_k"]=22.6628895184
    bad=case_screen(case)
    assert "wall_coefficient_inconsistent_with_case_geometry" in bad["open_engineering_evidence_gates"]


def test_ignored_legacy_schedule_timing_is_an_open_qualification_gate():
    current = case_screen(wang_case())
    assert current["schedule_timing_inputs_used"] == {"slot_s": 120.0}
    assert current["schedule_timing_inputs_ignored"] == {}
    legacy = wang_case()
    legacy["schedule"]["pressurization_s"] = 30.0
    screened = case_screen(legacy)
    assert screened["schedule_timing_inputs_ignored"] == {"pressurization_s": 30.0}
    assert "schedule_timing_inputs_ignored_by_plugin" in screened["open_engineering_evidence_gates"]


def test_conflicting_declared_cycle_time_cannot_pass_input_screen():
    case = wang_case()
    case["schedule"]["declared_cycle_time_s"] = 1500.0
    screened = case_screen(case)
    assert "declared_cycle_time_conflicts_with_executed_schedule" in screened["open_engineering_evidence_gates"]
