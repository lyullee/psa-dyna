"""Convert a CSS cycle ledger into a guarded dynamic-PSA mean stream.

The bridge is intentionally stricter than a normal post-processor.  It only
converts cycle-integrated feed/product/tail ledgers after an ordinary CSS run,
and it never authorizes DWSIM, TEA/LCA, design, or training use unless the
upstream source and material-domain gates are explicitly open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .release_filter import evaluate as evaluate_release
from .phase_aligned_state import resolve_phase_aligned_states
from .schedules import build_schedule


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _finite_nonnegative(values: list[Any], label: str) -> list[float]:
    output = [float(value) for value in values]
    if not output or any(value < -1e-12 for value in output):
        raise ValueError(f"{label} must be finite and non-negative")
    return [max(value, 0.0) for value in output]


def _mole_fractions(rates: list[float]) -> list[float]:
    """Return a finite composition without inventing a composition for zero flow."""
    total = float(sum(rates))
    if total <= 0.0:
        raise ValueError("A mean-stream composition requires positive total molar flow")
    return [float(value / total) for value in rates]


def _property_reference_complete(reference: Any) -> bool:
    if not isinstance(reference, dict):
        return False
    if any(
        not isinstance(reference.get(name), str) or not reference[name].strip()
        for name in ("name", "version", "source", "sha256")
    ):
        return False
    checksum = reference["sha256"].strip().lower()
    return len(checksum) == 64 and all(character in "0123456789abcdef" for character in checksum)


def _numerical_resolution_evidence(
    result_path: Path,
    result: dict[str, Any],
    audit_paths: dict[str, Path] | None,
) -> dict[str, Any]:
    """Require hash-linked space, time and dynamic-profile qualification.

    A CSS result is not automatically mesh independent.  The bridge therefore
    treats the three numerical audits as explicit upstream evidence instead of
    inferring qualification from a case name or cell count.
    """
    required = ("spatial", "temporal", "profile")
    if not audit_paths or any(name not in audit_paths for name in required):
        return {
            "status": "missing",
            "qualified": False,
            "errors": ["spatial_temporal_profile_audits_required"],
            "audits": {},
        }
    result_hash = digest(result_path).lower()
    engine_hash = str(result.get("code_sha256") or "").lower()
    records: dict[str, Any] = {}
    errors: list[str] = []
    payloads: dict[str, dict[str, Any]] = {}
    for name in required:
        path = Path(audit_paths[name]).resolve()
        if not path.is_file():
            errors.append(f"{name}_audit_missing")
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            errors.append(f"{name}_audit_unreadable")
            continue
        payloads[name] = payload
        records[name] = {"path": str(path), "sha256": digest(path), "status": payload.get("status")}

    spatial = payloads.get("spatial", {})
    spatial_hashes = {
        str(row.get("sha256") or "").lower()
        for row in spatial.get("inputs", [])
        if isinstance(row, dict)
    }
    spatial_pass = bool(
        spatial.get("status") == "three_grid_screen_pass_not_source_release"
        and spatial.get("numerical_screen_pass") is True
        and result_hash in spatial_hashes
    )
    if not spatial_pass:
        errors.append("spatial_audit_failed_or_result_unlinked")

    temporal = payloads.get("temporal", {})
    temporal_pass = bool(
        temporal.get("status") == "diagnostic_time_refinement_pass_not_source_release"
        and (temporal.get("pairwise_numerical_screen") or {}).get("pass") is True
        and temporal.get("same_engine_hash") is True
        and engine_hash
        and str(temporal.get("engine_code_sha256") or "").lower() == engine_hash
    )
    if not temporal_pass:
        errors.append("temporal_audit_failed_or_engine_unlinked")

    profile = payloads.get("profile", {})
    profile_hashes = {
        str(row.get("result_sha256") or "").lower()
        for row in profile.get("inputs", [])
        if isinstance(row, dict)
    }
    profile_pass = bool(
        profile.get("status") == "profile_grid_screen_pass_not_source_release"
        and profile.get("numerical_screen_pass") is True
        and result_hash in profile_hashes
    )
    if not profile_pass:
        errors.append("profile_audit_failed_or_result_unlinked")

    qualified = not errors
    return {
        "status": "qualified" if qualified else "failed_closed",
        "qualified": qualified,
        "result_sha256": result_hash,
        "engine_code_sha256": engine_hash or None,
        "errors": sorted(set(errors)),
        "audits": records,
        "interpretation": (
            "All three numerical audits pass and the spatial/profile evidence names this exact result hash; "
            "the temporal evidence names the same engine hash. This does not open source, design or DWSIM gates."
            if qualified
            else "CSS alone is insufficient; missing, failed or unlinked numerical evidence keeps downstream use closed."
        ),
    }


def build_mean_stream(
    result_path: Path,
    case_path: Path,
    output_path: Path | None = None,
    qualification_audits: dict[str, Path] | None = None,
) -> dict[str, Any]:
    result_path = Path(result_path).resolve()
    case_path = Path(case_path).resolve()
    result = json.loads(result_path.read_text(encoding="utf-8-sig"))
    case = json.loads(case_path.read_text(encoding="utf-8-sig"))
    case_hash = digest(case_path)
    # Use the common fail-closed release filter rather than relying only on
    # legacy fields such as dataset_release_gate.  This keeps the bridge safe
    # if a future result schema changes its gate names or adds contradictory
    # flags.  The decisions are recorded in the output for auditability.
    upstream_release_filter = {
        purpose: evaluate_release(result_path, purpose)
        for purpose in ("source_dynamic", "design", "dwsim", "tea_lca")
    }
    numerical_resolution = _numerical_resolution_evidence(
        result_path, result, qualification_audits
    )
    phase_aligned_state = resolve_phase_aligned_states(result)
    if result.get("input_sha256", "").lower() != case_hash.lower():
        raise ValueError("Result input_sha256 does not match the supplied case")

    final = result.get("final_cycle") or {}
    status = result.get("status")
    css_streak = int(final.get("css_streak", 0) or 0)
    state_css = bool(result.get("state_css_converged", False))
    solver = case.get("solver", {})
    if status != "css_converged" or not state_css or css_streak < int(solver.get("css_consecutive_cycles", 3)):
        raise ValueError("Only an ordinary css_converged result with the required CSS streak may be transferred")

    mass_relative = float(final.get("mass_balance_relative"))
    energy_relative = float(final.get("energy_balance_relative"))
    if mass_relative > float(solver.get("mass_balance_tolerance", 1e-5)):
        raise ValueError("Mass-balance residual exceeds the case tolerance")
    if energy_relative > float(solver.get("energy_balance_tolerance", 1e-5)):
        raise ValueError("Energy-balance residual exceeds the case tolerance")

    schedule = build_schedule(case["schedule"])
    cycle_time_s = float(schedule["cycle_time_s"])
    if cycle_time_s <= 0.0:
        raise ValueError("Cycle time must be positive")
    species = [item["name"] for item in case["species"]]
    feed = _finite_nonnegative(final.get("feed_by_component_mol", []), "feed_by_component_mol")
    product = _finite_nonnegative(final.get("product_by_component_mol", []), "product_by_component_mol")
    tail = _finite_nonnegative(final.get("tail_by_component_mol", []), "tail_by_component_mol")
    if not (len(feed) == len(product) == len(tail) == len(species)):
        raise ValueError("Feed, product, tail and species lengths do not match")

    feed_rate = [value / cycle_time_s for value in feed]
    product_rate = [value / cycle_time_s for value in product]
    tail_rate = [value / cycle_time_s for value in tail]
    h2_index = species.index("H2")
    product_total = sum(product_rate)
    feed_h2 = feed_rate[h2_index]
    h2_purity = product_rate[h2_index] / product_total if product_total > 0.0 else None
    h2_recovery = product_rate[h2_index] / feed_h2 if feed_h2 > 0.0 else None
    if h2_purity is None or h2_recovery is None:
        raise ValueError("H2 product and feed rates must be positive")
    if abs(h2_purity - float(final["h2_purity"])) > 1e-8:
        raise ValueError("Ledger-derived H2 purity is inconsistent with the result record")
    if abs(h2_recovery - float(final["h2_recovery"])) > 1e-8:
        raise ValueError("Ledger-derived H2 recovery is inconsistent with the result record")

    # New result records carry the three boundary enthalpy ledgers and the
    # vacuum-work ledger explicitly.  Older archived results predate this
    # schema extension; they remain readable as numerical mean streams, but
    # they are deliberately not eligible for an enthalpy-aware DWSIM bridge.
    enthalpy_keys = ("feed_enthalpy_j", "product_enthalpy_j", "tail_enthalpy_j", "ambient_heat_j", "vacuum_work_j")
    enthalpy_available = all(key in final for key in enthalpy_keys)
    enthalpy_ledger = None
    if enthalpy_available:
        enthalpy_ledger = {}
        for key in enthalpy_keys:
            value = float(final[key])
            if not (value == value and abs(value) != float("inf")):
                raise ValueError(f"{key} must be finite")
            enthalpy_ledger[key] = value
        if enthalpy_ledger["vacuum_work_j"] < -1e-9:
            raise ValueError("vacuum_work_j must be non-negative")
        enthalpy_ledger["energy_balance_relative"] = float(final.get("energy_balance_relative", float("nan")))
        if not (enthalpy_ledger["energy_balance_relative"] == enthalpy_ledger["energy_balance_relative"]):
            raise ValueError("energy_balance_relative must be finite when enthalpy ledger is present")
        enthalpy_ledger["feed_enthalpy_w"] = enthalpy_ledger["feed_enthalpy_j"] / cycle_time_s
        enthalpy_ledger["product_enthalpy_w"] = enthalpy_ledger["product_enthalpy_j"] / cycle_time_s
        enthalpy_ledger["tail_enthalpy_w"] = enthalpy_ledger["tail_enthalpy_j"] / cycle_time_s
        enthalpy_ledger["ambient_heat_w"] = enthalpy_ledger["ambient_heat_j"] / cycle_time_s
        enthalpy_ledger["vacuum_work_w"] = enthalpy_ledger["vacuum_work_j"] / cycle_time_s
        enthalpy_ledger["feed_specific_enthalpy_j_mol"] = enthalpy_ledger["feed_enthalpy_j"] / max(sum(feed), 1e-30)
        enthalpy_ledger["product_specific_enthalpy_j_mol"] = enthalpy_ledger["product_enthalpy_j"] / max(sum(product), 1e-30)
        enthalpy_ledger["tail_specific_enthalpy_j_mol"] = enthalpy_ledger["tail_enthalpy_j"] / max(sum(tail), 1e-30)

    energy_mode = case.get("bed", {}).get("energy_mode", "unspecified")
    enthalpy_basis = {
        "status": "available" if enthalpy_available else "missing_in_archived_result",
        "energy_mode": energy_mode,
        "definition": final.get("energy_ledger_basis", "not recorded"),
        "property_package": "not assigned; model ledger only",
        "transfer_ready": False,
    }
    if enthalpy_ledger is not None:
        enthalpy_basis["ledger"] = enthalpy_ledger

    # DWSIM requires a thermodynamic state in addition to a composition and
    # flow.  A cycle-average ledger does not define one scalar pressure or
    # temperature, so the payload deliberately leaves those values null.  A
    # downstream adapter must obtain them from a phase-aligned ledger or a
    # DWSIM/property-package flash; this prevents accidental use of the case
    # boundary pressure as a fictitious mean stream pressure.
    enthalpy_flow = {
        "FEED": (None if enthalpy_ledger is None else enthalpy_ledger["feed_enthalpy_w"]),
        "PRODUCT": (None if enthalpy_ledger is None else enthalpy_ledger["product_enthalpy_w"]),
        "TAIL": (None if enthalpy_ledger is None else enthalpy_ledger["tail_enthalpy_w"]),
    }
    dwsim_streams = {}
    for name, rates in (("FEED", feed_rate), ("PRODUCT", product_rate), ("TAIL", tail_rate)):
        state = phase_aligned_state["streams"][name]
        state_resolved = state.get("status") == "resolved"
        dwsim_streams[name] = {
            "molar_flow_mol_s": float(sum(rates)),
            "mole_fractions": _mole_fractions(rates),
            "pressure_bar": state.get("pressure_bar") if state_resolved else None,
            "temperature_k": state.get("temperature_k") if state_resolved else None,
            "enthalpy_flow_w": enthalpy_flow[name],
            "state_basis": (state.get("basis") if state_resolved else
                            "cycle_average_composition_and_ledger; phase-aligned scalar P/T unresolved"),
        }

    trajectory_in_domain = result.get("trajectory_domain_sampled_within_documented") is True
    source_gate_open = upstream_release_filter["source_dynamic"]["accepted"]
    design_gate_open = upstream_release_filter["design"]["accepted"]
    dwsim_gate_open = upstream_release_filter["dwsim"]["accepted"]
    tea_lca_gate_open = upstream_release_filter["tea_lca"]["accepted"]
    # A source/material pass is necessary but not sufficient for a DWSIM
    # transfer.  Scalar P/T are accepted only when all three streams contain
    # explicit phase-aligned samples; otherwise the bridge remains closed.
    dwsim_state_resolved = phase_aligned_state["status"] == "resolved"
    property_package_reference = case.get("dwsim_property_package_reference")
    if not isinstance(property_package_reference, dict):
        property_package_reference = None
    property_reference_complete = _property_reference_complete(property_package_reference)
    common_transfer_ready = bool(
        trajectory_in_domain
        and numerical_resolution["qualified"]
        and source_gate_open
        and design_gate_open
        and enthalpy_available
        and dwsim_state_resolved
    )
    dwsim_transfer_allowed = bool(
        common_transfer_ready and dwsim_gate_open and property_reference_complete
    )
    tea_lca_transfer_allowed = bool(common_transfer_ready and tea_lca_gate_open)
    dwsim_blockers: list[str] = []
    if not dwsim_state_resolved:
        dwsim_blockers.append("phase_aligned_scalar_pressure_temperature_unresolved")
    if not numerical_resolution["qualified"]:
        dwsim_blockers.append(
            "numerical_spatial_temporal_profile_qualification_missing_or_failed"
        )
    if not property_reference_complete:
        dwsim_blockers.append("validated_dwsim_property_package_unassigned_or_incomplete")
    if not source_gate_open:
        dwsim_blockers.append("source_dynamic_validation_not_released")
    if not design_gate_open:
        dwsim_blockers.append("design_qualification_not_released")
    if not dwsim_gate_open:
        dwsim_blockers.append("upstream_dwsim_release_not_authorized")
    report = {
        "schema": "dynamic-psa-mean-stream/v1",
        "status": "qualified_numerical_mean_stream_closed_for_downstream_release",
        "result": {"path": str(result_path), "sha256": digest(result_path)},
        "case": {"path": str(case_path), "sha256": case_hash, "case_id": case.get("case_id")},
        "topology": {
            "number_of_beds": result.get("number_of_beds"),
            "number_of_cycle_steps": result.get("number_of_cycle_steps"),
            "axial_cells": case.get("bed", {}).get("n_cells"),
        },
        "css_qualification": {
            "status": status,
            "state_css_converged": state_css,
            "css_streak": css_streak,
            "state_scaled_max": final.get("state_scaled_max"),
            "mass_balance_relative": mass_relative,
            "energy_balance_relative": energy_relative,
        },
        "mean_stream_basis": {
            "definition": "cycle-integrated ordinary-run ledgers divided by schedule cycle_time_s",
            "cycle_time_s": cycle_time_s,
            "species_order": species,
            "feed_mol_s": feed_rate,
            "product_mol_s": product_rate,
            "tail_mol_s": tail_rate,
            "h2_purity": h2_purity,
            "h2_recovery": h2_recovery,
            "vacuum_work_w": float(final.get("vacuum_work_j", 0.0)) / cycle_time_s,
        },
        "enthalpy_basis": enthalpy_basis,
        "dwsim_input_payload": {
            "schema": "dynamic-psa-dwsim-mean-stream-input/v1",
            "status": (
                "validated_ready_for_dwsim"
                if dwsim_transfer_allowed
                else "prepared_closed_pending_qualification"
            ),
            "components": species,
            "streams": dwsim_streams,
            "units": {
                "molar_flow": "mol/s",
                "mole_fraction": "1",
                "pressure": "bar",
                "temperature": "K",
                "enthalpy_flow": "W",
            },
            "pressure_temperature_policy": "Do not infer scalar cycle-average pressure or temperature from boundary values; resolve with phase-aligned stream data and a validated DWSIM property package.",
            "property_package": property_package_reference,
            "property_package_reference_complete": property_reference_complete,
            "dwsim_state_resolved": dwsim_state_resolved,
            "numerical_resolution_qualified": numerical_resolution["qualified"],
            "release_blockers": dwsim_blockers,
            "release_flags": {
                "dwsim_transfer_allowed": dwsim_transfer_allowed,
                "tea_lca_transfer_allowed": tea_lca_transfer_allowed,
            },
        },
        "release_flags": {
            "numerical_mean_stream": True,
            "enthalpy_aware_mean_stream": enthalpy_available,
            "dwsim_transfer_allowed": dwsim_transfer_allowed,
            "tea_lca_transfer_allowed": tea_lca_transfer_allowed,
            "design_eligible": bool(common_transfer_ready),
            "training_eligible": bool(common_transfer_ready and source_gate_open),
        },
        "upstream_gate_observations": {
            "trajectory_domain_sampled_within_documented": result.get("trajectory_domain_sampled_within_documented"),
            "dataset_release_gate": result.get("dataset_release_gate"),
            "source_code_zip_sha256": result.get("source_code_zip_sha256"),
        },
        "upstream_release_filter": upstream_release_filter,
        "numerical_resolution_qualification": numerical_resolution,
        "phase_aligned_state": phase_aligned_state,
        "claim_boundary": "This is a numerically qualified cycle-average stream. The enthalpy ledger is a model energy ledger, not a validated DWSIM property-package enthalpy. It is not source validation or a plant-design stream until the explicit upstream source/material-domain and property-reference gates are open.",
    }
    if output_path is not None:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spatial-audit", type=Path)
    parser.add_argument("--temporal-audit", type=Path)
    parser.add_argument("--profile-audit", type=Path)
    args = parser.parse_args()
    audit_values = {
        "spatial": args.spatial_audit,
        "temporal": args.temporal_audit,
        "profile": args.profile_audit,
    }
    if any(audit_values.values()) and not all(audit_values.values()):
        parser.error(
            "--spatial-audit, --temporal-audit and --profile-audit must be supplied together"
        )
    qualification_audits = (
        audit_values if all(audit_values.values()) else None
    )
    report = build_mean_stream(
        args.result,
        args.case,
        args.output,
        qualification_audits=qualification_audits,
    )
    print(json.dumps({"schema": report["schema"], "status": report["status"], "output": str(args.output), "release_flags": report["release_flags"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
