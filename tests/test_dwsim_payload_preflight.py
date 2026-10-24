import json
from pathlib import Path

from psa_dynamic.dwsim_payload_preflight import validate_payload


def test_structurally_valid_payload_remains_closed_for_dwsim(tmp_path):
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2", "CO2"],
        "streams": {
            name: {
                "molar_flow_mol_s": 1.0,
                "mole_fractions": [0.9, 0.1],
                "pressure_bar": None,
                "temperature_k": None,
                "enthalpy_flow_w": 10.0,
            }
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "dwsim_state_resolved": False,
        "numerical_resolution_qualified": True,
        "status": "prepared_closed_pending_qualification",
    }
    path = tmp_path / "test_dwsim_input_payload.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert report["status"] == "contract_valid_closed_pending_property_reference"
    assert report["integrity_checks"]["structural_pass"] is True
    assert len(report["payload"]["sha256"]) == 64
    assert report["release_flags"]["dwsim_transfer_allowed"] is False
    assert "FEED:pressure_unresolved" in report["blockers"]
    assert "property_package_unassigned" in report["blockers"]


def test_invalid_composition_is_rejected(tmp_path):
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2", "CO2"],
        "streams": {
            name: {
                "molar_flow_mol_s": 1.0,
                "mole_fractions": [0.8, 0.3],
                "pressure_bar": None,
                "temperature_k": None,
                "enthalpy_flow_w": 1.0,
            }
            for name in ("FEED", "PRODUCT", "TAIL")
        },
    }
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert report["status"] == "contract_rejected"
    assert "FEED:mole_fractions_do_not_sum_to_one" in report["errors"]


def test_property_package_requires_hash_pinned_reference(tmp_path):
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2"],
        "streams": {
            name: {
                "molar_flow_mol_s": 1.0,
                "mole_fractions": [1.0],
                "pressure_bar": 10.0,
                "temperature_k": 300.0,
                "enthalpy_flow_w": 1.0,
            }
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "property_package": {
            "name": "Peng-Robinson",
            "version": "example-1",
            "source": "local-test",
            "sha256": "0" * 64,
        },
        "dwsim_state_resolved": True,
        "numerical_resolution_qualified": True,
        "status": "validated_ready_for_dwsim",
    }
    path = tmp_path / "complete.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert report["status"] == "contract_valid_closed_pending_property_reference"
    assert "property_package_reference_incomplete" not in report["blockers"]
    assert "property_package_name_missing" not in report["blockers"]
    assert "numerical_resolution_not_qualified" not in report["blockers"]


def test_property_package_string_is_not_a_qualified_reference(tmp_path):
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2"],
        "streams": {
            name: {"molar_flow_mol_s": 1.0, "mole_fractions": [1.0],
                   "pressure_bar": 10.0, "temperature_k": 300.0}
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "property_package": "Peng-Robinson",
    }
    path = tmp_path / "legacy.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert "property_package_reference_incomplete" in report["blockers"]


def test_property_package_file_hash_is_verified_when_path_is_supplied(tmp_path):
    package_file = tmp_path / "property_package.dll"
    package_file.write_bytes(b"qualified-test-package")
    import hashlib
    package_hash = hashlib.sha256(package_file.read_bytes()).hexdigest()
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2"],
        "streams": {
            name: {"molar_flow_mol_s": 1.0, "mole_fractions": [1.0],
                   "pressure_bar": 10.0, "temperature_k": 300.0}
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "property_package": {
            "name": "test-package", "version": "1", "source": "fixture",
            "path": package_file.name, "sha256": package_hash,
        },
    }
    path = tmp_path / "payload.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert report["property_package"]["file_hash_verified"] is True
    assert "property_package_sha256_mismatch" not in report["blockers"]


def test_property_package_file_hash_mismatch_is_blocked(tmp_path):
    package_file = tmp_path / "property_package.dll"
    package_file.write_bytes(b"actual-package")
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2"],
        "streams": {
            name: {"molar_flow_mol_s": 1.0, "mole_fractions": [1.0],
                   "pressure_bar": 10.0, "temperature_k": 300.0}
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "property_package": {
            "name": "test-package", "version": "1", "source": "fixture",
            "path": package_file.name, "sha256": "0" * 64,
        },
    }
    path = tmp_path / "payload.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert "property_package_sha256_mismatch" in report["blockers"]


def test_prepared_effective_closure_cannot_be_misread_as_transfer_release(tmp_path):
    package_file = tmp_path / "OpenPRMHV2.DWSIM.dll"
    package_file.write_bytes(b"effective-closure-fixture")
    import hashlib
    payload = {
        "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
        "components": ["H2"],
        "streams": {
            name: {
                "molar_flow_mol_s": 1.0,
                "mole_fractions": [1.0],
                "pressure_bar": 10.0,
                "temperature_k": 300.0,
            }
            for name in ("FEED", "PRODUCT", "TAIL")
        },
        "property_package": {
            "name": "Peng-Robinson-MHV2 (Open Prototype)",
            "version": "DWSIM 8.7.1 fixture",
            "source": "effective fitted closure; not native Aspen equivalence",
            "path": package_file.name,
            "sha256": hashlib.sha256(package_file.read_bytes()).hexdigest(),
            "qualification": {"eligible_for_payload_injection": False},
        },
        "dwsim_state_resolved": True,
        "numerical_resolution_qualified": True,
        "status": "validated_ready_for_dwsim",
    }
    path = tmp_path / "payload.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    report = validate_payload(path)
    assert report["property_package"]["file_hash_verified"] is True
    assert report["property_package"]["eligible_for_payload_injection"] is False
    assert "property_package_reference_not_eligible_for_payload_injection" in report["blockers"]
    assert report["release_flags"]["dwsim_transfer_allowed"] is False
