from __future__ import annotations

import hashlib
import json

import pytest

from psa_dynamic.mean_stream_bridge import build_mean_stream


def _write_case_and_result(tmp_path, status="css_converged"):
    case_path = tmp_path / "case.json"
    case = {
        "case_id": "bridge_case",
        "species": [{"name": "H2"}, {"name": "CH4"}],
        "schedule": {"name": "skarstrom_4step", "pressurization_s": 1.0, "adsorption_s": 1.0},
        "solver": {
            "css_consecutive_cycles": 3,
            "mass_balance_tolerance": 1e-5,
            "energy_balance_tolerance": 1e-5,
        },
        "bed": {"n_cells": 2, "energy_mode": "internal_energy"},
    }
    case_path.write_text(json.dumps(case), encoding="utf-8")
    digest = hashlib.sha256(case_path.read_bytes()).hexdigest()
    result_path = tmp_path / "result.json"
    result = {
        "schema": "dynamic-psa-result/v1",
        "status": status,
        "input_sha256": digest,
        "code_sha256": "engine-hash",
        "state_css_converged": status == "css_converged",
        "number_of_beds": 2,
        "number_of_cycle_steps": 4,
        "dataset_release_gate": "closed_pending_independent_resolution_and_source_validation_audit",
        "trajectory_domain_sampled_within_documented": False,
        "final_cycle": {
            "css_streak": 3 if status == "css_converged" else 0,
            "state_scaled_max": 1e-5,
            "mass_balance_relative": 1e-8,
            "energy_balance_relative": 1e-8,
            "feed_by_component_mol": [4.0, 2.0],
            "product_by_component_mol": [3.0, 0.5],
            "tail_by_component_mol": [1.0, 1.5],
            "h2_purity": 3.0 / 3.5,
            "h2_recovery": 3.0 / 4.0,
            "feed_enthalpy_j": 100.0,
            "product_enthalpy_j": 70.0,
            "tail_enthalpy_j": 20.0,
            "ambient_heat_j": 2.0,
            "vacuum_work_j": 8.0,
            "energy_ledger_basis": "internal_energy: sensible enthalpy ledger plus vacuum work",
        },
    }
    result_path.write_text(json.dumps(result), encoding="utf-8")
    return case_path, result_path


def test_bridge_converts_cycle_ledgers_and_keeps_downstream_closed(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    report = build_mean_stream(result_path, case_path)
    assert report["mean_stream_basis"]["cycle_time_s"] == 4.0
    assert report["mean_stream_basis"]["product_mol_s"] == [0.75, 0.125]
    assert report["mean_stream_basis"]["h2_purity"] == pytest.approx(3.0 / 3.5)
    assert report["mean_stream_basis"]["h2_recovery"] == pytest.approx(0.75)
    assert report["release_flags"]["numerical_mean_stream"] is True
    assert report["release_flags"]["enthalpy_aware_mean_stream"] is True
    assert report["enthalpy_basis"]["status"] == "available"
    assert report["enthalpy_basis"]["ledger"]["product_enthalpy_w"] == pytest.approx(17.5)
    assert report["enthalpy_basis"]["ledger"]["product_specific_enthalpy_j_mol"] == pytest.approx(20.0)
    assert report["release_flags"]["dwsim_transfer_allowed"] is False
    assert report["release_flags"]["tea_lca_transfer_allowed"] is False
    assert report["upstream_release_filter"]["source_dynamic"]["accepted"] is False
    assert report["upstream_release_filter"]["dwsim"]["accepted"] is False
    assert report["numerical_resolution_qualification"]["status"] == "missing"
    assert report["dwsim_input_payload"]["numerical_resolution_qualified"] is False
    assert "numerical_spatial_temporal_profile_qualification_missing_or_failed" in report["dwsim_input_payload"]["release_blockers"]


def test_bridge_rejects_non_css_result(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path, status="maximum_cycles_reached")
    with pytest.raises(ValueError, match="css_converged"):
        build_mean_stream(result_path, case_path)


def test_bridge_keeps_archived_result_closed_when_enthalpy_ledger_is_missing(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["final_cycle"].pop("feed_enthalpy_j")
    result["final_cycle"].pop("product_enthalpy_j")
    result["final_cycle"].pop("tail_enthalpy_j")
    result_path.write_text(json.dumps(result), encoding="utf-8")
    report = build_mean_stream(result_path, case_path)
    assert report["enthalpy_basis"]["status"] == "missing_in_archived_result"
    assert report["release_flags"]["enthalpy_aware_mean_stream"] is False
    assert report["release_flags"]["dwsim_transfer_allowed"] is False


def test_bridge_does_not_treat_legacy_dataset_gate_as_release_permission(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["dataset_release_gate"] = "open"
    result_path.write_text(json.dumps(result), encoding="utf-8")
    report = build_mean_stream(result_path, case_path)
    assert report["upstream_release_filter"]["source_dynamic"]["accepted"] is False
    assert report["release_flags"]["dwsim_transfer_allowed"] is False


def test_bridge_uses_phase_aligned_states_when_explicitly_present(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    samples = [
        {"time_s": 0.0, "duration_s": 1.0, "flow_mol_s": 1.0,
         "pressure_bar": 8.0, "temperature_k": 300.0},
    ]
    result["phase_aligned_streams"] = {
        "FEED": samples, "PRODUCT": samples, "TAIL": samples,
    }
    result_path.write_text(json.dumps(result), encoding="utf-8")
    report = build_mean_stream(result_path, case_path)
    assert report["phase_aligned_state"]["status"] == "resolved"
    assert report["dwsim_input_payload"]["dwsim_state_resolved"] is True
    assert report["dwsim_input_payload"]["streams"]["PRODUCT"]["pressure_bar"] == 8.0
    assert report["release_flags"]["dwsim_transfer_allowed"] is False
    assert "phase_aligned_scalar_pressure_temperature_unresolved" not in report["dwsim_input_payload"]["release_blockers"]


def test_bridge_propagates_hash_pinned_property_package_reference(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    case = json.loads(case_path.read_text(encoding="utf-8"))
    case["dwsim_property_package_reference"] = {
        "name": "Peng-Robinson",
        "version": "DWSIM-test",
        "source": "local-test-reference",
        "sha256": "1" * 64,
    }
    case_path.write_text(json.dumps(case), encoding="utf-8")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["input_sha256"] = hashlib.sha256(case_path.read_bytes()).hexdigest()
    result_path.write_text(json.dumps(result), encoding="utf-8")
    report = build_mean_stream(result_path, case_path)
    assert report["dwsim_input_payload"]["property_package"]["name"] == "Peng-Robinson"
    assert report["dwsim_input_payload"]["property_package_reference_complete"] is True
    assert "validated_dwsim_property_package_unassigned_or_incomplete" not in report["dwsim_input_payload"]["release_blockers"]
    assert report["release_flags"]["dwsim_transfer_allowed"] is False


def _write_resolution_audits(tmp_path, result_path, *, linked=True):
    result_hash = hashlib.sha256(result_path.read_bytes()).hexdigest()
    named_hash = result_hash if linked else "0" * 64
    spatial = tmp_path / "spatial.json"
    temporal = tmp_path / "temporal.json"
    profile = tmp_path / "profile.json"
    spatial.write_text(json.dumps({
        "status": "three_grid_screen_pass_not_source_release",
        "numerical_screen_pass": True,
        "inputs": [{"sha256": named_hash}],
    }), encoding="utf-8")
    temporal.write_text(json.dumps({
        "status": "diagnostic_time_refinement_pass_not_source_release",
        "pairwise_numerical_screen": {"pass": True},
        "same_engine_hash": True,
        "engine_code_sha256": "engine-hash",
    }), encoding="utf-8")
    profile.write_text(json.dumps({
        "status": "profile_grid_screen_pass_not_source_release",
        "numerical_screen_pass": True,
        "inputs": [{"result_sha256": named_hash}],
    }), encoding="utf-8")
    return {"spatial": spatial, "temporal": temporal, "profile": profile}


def test_bridge_records_hash_linked_numerical_resolution_evidence(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    audits = _write_resolution_audits(tmp_path, result_path)
    report = build_mean_stream(
        result_path, case_path, qualification_audits=audits
    )
    assert report["numerical_resolution_qualification"]["qualified"] is True
    assert report["dwsim_input_payload"]["numerical_resolution_qualified"] is True
    assert "numerical_spatial_temporal_profile_qualification_missing_or_failed" not in report["dwsim_input_payload"]["release_blockers"]
    # Numerical qualification is necessary but still not sufficient: source,
    # design, phase-aligned state and property-package gates remain separate.
    assert report["release_flags"]["dwsim_transfer_allowed"] is False


def test_bridge_rejects_resolution_audits_not_linked_to_exact_result(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    audits = _write_resolution_audits(tmp_path, result_path, linked=False)
    report = build_mean_stream(
        result_path, case_path, qualification_audits=audits
    )
    qualification = report["numerical_resolution_qualification"]
    assert qualification["qualified"] is False
    assert "spatial_audit_failed_or_result_unlinked" in qualification["errors"]
    assert "profile_audit_failed_or_result_unlinked" in qualification["errors"]


def test_bridge_has_a_real_future_ready_path_when_every_independent_gate_is_open(tmp_path):
    case_path, result_path = _write_case_and_result(tmp_path)
    case = json.loads(case_path.read_text(encoding="utf-8"))
    case["dwsim_property_package_reference"] = {
        "name": "Peng-Robinson",
        "version": "qualified-test-version",
        "source": "hash-pinned-test-reference",
        "sha256": "a" * 64,
    }
    case_path.write_text(json.dumps(case), encoding="utf-8")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["input_sha256"] = hashlib.sha256(case_path.read_bytes()).hexdigest()
    result["trajectory_domain_sampled_within_documented"] = True
    result["release_flags"] = {
        "source_dynamic_validation": True,
        "design_eligible": True,
        "dwsim_application_allowed": True,
        "tea_lca_allowed": True,
    }
    samples = [{
        "time_s": 0.0,
        "duration_s": 1.0,
        "flow_mol_s": 1.0,
        "pressure_bar": 8.0,
        "temperature_k": 300.0,
    }]
    result["phase_aligned_streams"] = {
        "FEED": samples,
        "PRODUCT": samples,
        "TAIL": samples,
    }
    result_path.write_text(json.dumps(result), encoding="utf-8")
    audits = _write_resolution_audits(tmp_path, result_path)
    report = build_mean_stream(
        result_path, case_path, qualification_audits=audits
    )
    assert report["release_flags"]["dwsim_transfer_allowed"] is True
    assert report["dwsim_input_payload"]["status"] == "validated_ready_for_dwsim"
    assert report["dwsim_input_payload"]["release_blockers"] == []
