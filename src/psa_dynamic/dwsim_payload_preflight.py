"""Read-only structural preflight for a dynamic-PSA DWSIM payload.

This module checks the handoff contract without opening DWSIM and without
authorizing a flowsheet transfer.  A cycle-average stream is not a complete
thermodynamic state until phase-aligned pressure, temperature and a named
property package have been supplied.  The default result therefore remains a
closed preparation result even when the JSON is structurally valid.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


EXPECTED_SCHEMA = "dynamic-psa-dwsim-mean-stream-input/v1"
EXPECTED_STREAMS = ("FEED", "PRODUCT", "TAIL")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _fractions(values: Any, count: int, label: str, errors: list[str]) -> None:
    if not isinstance(values, list) or len(values) != count:
        errors.append(f"{label}:mole_fractions_length")
        return
    converted: list[float] = []
    for index, value in enumerate(values):
        if not _finite(value) or float(value) < 0.0:
            errors.append(f"{label}:mole_fraction_{index}_invalid")
        else:
            converted.append(float(value))
    if converted and abs(sum(converted) - 1.0) > 1e-8:
        errors.append(f"{label}:mole_fractions_do_not_sum_to_one")


def validate_payload(path: Path) -> dict[str, Any]:
    """Validate one bridge payload and return a JSON-serializable audit."""
    path = Path(path).resolve()
    errors: list[str] = []
    blockers: list[str] = []
    try:
        raw_payload = json.loads(path.read_text(encoding="utf-8-sig"))
        # Accept either the standalone dwsim_input_payload object or the
        # complete mean-stream bridge that contains it.  This lets a user
        # preflight the exact archived bridge without manually extracting a
        # nested JSON object.
        payload = raw_payload.get("dwsim_input_payload", raw_payload)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return {
            "schema": "dynamic-psa-dwsim-payload-preflight/v1",
            "status": "contract_rejected",
            "payload": {"path": str(path), "exists": path.is_file(), "sha256": sha256(path) if path.is_file() else None},
            "errors": [f"payload_read_error:{type(error).__name__}"],
            "blockers": [],
            "release_flags": {"dwsim_transfer_allowed": False},
        }

    if payload.get("schema") != EXPECTED_SCHEMA:
        errors.append("schema_mismatch")
    components = payload.get("components")
    if not isinstance(components, list) or not components:
        errors.append("components_missing_or_empty")
        components = []
    elif any(not isinstance(name, str) or not name for name in components):
        errors.append("component_name_invalid")
    elif len(set(components)) != len(components):
        errors.append("component_names_not_unique")

    streams = payload.get("streams")
    if not isinstance(streams, dict):
        errors.append("streams_missing_or_invalid")
        streams = {}
    else:
        for name in EXPECTED_STREAMS:
            if name not in streams:
                errors.append(f"stream_missing:{name}")
                continue
            stream = streams[name]
            if not isinstance(stream, dict):
                errors.append(f"stream_invalid:{name}")
                continue
            flow = stream.get("molar_flow_mol_s")
            if not _finite(flow) or float(flow) <= 0.0:
                errors.append(f"{name}:molar_flow_invalid")
            _fractions(stream.get("mole_fractions"), len(components), name, errors)
            enthalpy = stream.get("enthalpy_flow_w")
            if enthalpy is not None and not _finite(enthalpy):
                errors.append(f"{name}:enthalpy_flow_invalid")
            if stream.get("pressure_bar") is None:
                blockers.append(f"{name}:pressure_unresolved")
            elif not _finite(stream.get("pressure_bar")) or float(stream["pressure_bar"]) <= 0.0:
                errors.append(f"{name}:pressure_invalid")
            if stream.get("temperature_k") is None:
                blockers.append(f"{name}:temperature_unresolved")
            elif not _finite(stream.get("temperature_k")) or float(stream["temperature_k"]) <= 0.0:
                errors.append(f"{name}:temperature_invalid")

    property_package = payload.get("property_package")
    property_package_audit: dict[str, Any] = {
        "reference_present": isinstance(property_package, dict),
        "file_hash_verified": False,
        "resolved_path": None,
    }
    if property_package is None:
        blockers.extend(("property_package_unassigned", "property_package_reference_incomplete"))
    elif not isinstance(property_package, dict):
        blockers.append("property_package_reference_incomplete")
    else:
        for field in ("name", "version", "source", "sha256"):
            value = property_package.get(field)
            if not isinstance(value, str) or not value.strip():
                blockers.append(f"property_package_{field}_missing")
        package_hash = property_package.get("sha256")
        if isinstance(package_hash, str) and package_hash.strip():
            normalized = package_hash.strip().lower()
            if len(normalized) != 64 or any(char not in "0123456789abcdef" for char in normalized):
                blockers.append("property_package_sha256_invalid")
            package_path_value = property_package.get("path")
            if package_path_value is not None:
                if not isinstance(package_path_value, str) or not package_path_value.strip():
                    blockers.append("property_package_path_invalid")
                else:
                    package_path = Path(package_path_value.strip())
                    if not package_path.is_absolute():
                        package_path = path.parent / package_path
                    package_path = package_path.resolve()
                    property_package_audit["resolved_path"] = str(package_path)
                    if not package_path.is_file():
                        blockers.append("property_package_file_missing")
                    elif len(normalized) == 64 and sha256(package_path).lower() != normalized:
                        blockers.append("property_package_sha256_mismatch")
                    elif len(normalized) == 64:
                        property_package_audit["file_hash_verified"] = True
        qualification = property_package.get("qualification")
        if isinstance(qualification, dict):
            property_package_audit["qualification_present"] = True
            property_package_audit["eligible_for_payload_injection"] = bool(
                qualification.get("eligible_for_payload_injection") is True
            )
            if qualification.get("eligible_for_payload_injection") is not True:
                blockers.append("property_package_reference_not_eligible_for_payload_injection")
    if payload.get("dwsim_state_resolved") is not True:
        blockers.append("dwsim_state_not_resolved")
    if payload.get("numerical_resolution_qualified") is not True:
        blockers.append("numerical_resolution_not_qualified")
    if payload.get("status") != "validated_ready_for_dwsim":
        blockers.append("payload_not_marked_validated_ready")

    structural_pass = not errors
    status = (
        "contract_valid_closed_pending_property_reference"
        if structural_pass
        else "contract_rejected"
    )
    return {
        "schema": "dynamic-psa-dwsim-payload-preflight/v1",
        "status": status,
        "payload": {
            "path": str(path),
            "sha256": sha256(path),
            "schema": payload.get("schema"),
            "components": components,
            "stream_names": sorted(streams),
        },
        "integrity_checks": {
            "structural_pass": structural_pass,
            "stream_count": len(streams),
            "expected_streams_present": all(name in streams for name in EXPECTED_STREAMS),
        },
        "property_package": property_package_audit,
        "errors": sorted(set(errors)),
        "blockers": sorted(set(blockers)),
        "release_flags": {"dwsim_transfer_allowed": False},
        "claim_boundary": "Structural handoff check only. This audit never authorizes DWSIM, TEA/LCA, design or training use.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("payload", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = validate_payload(args.payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"schema": report["schema"], "status": report["status"], "output": str(args.output), "release_flags": report["release_flags"]}, ensure_ascii=False))
    return 0 if report["integrity_checks"]["structural_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
