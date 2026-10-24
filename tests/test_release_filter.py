import json
from pathlib import Path

from psa_dynamic.release_filter import evaluate, filter_results


def write(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_default_deny_missing_and_false_flags(tmp_path):
    missing = tmp_path / "missing.json"
    false = write(tmp_path / "false.json", {"status": "css_converged"})
    report = filter_results([missing, false], "training")
    assert report["accepted_count"] == 0
    assert "file_missing" in report["rows"][0]["reasons"]
    assert "required_release_flag_missing" in report["rows"][1]["reasons"]


def test_true_flag_is_required_and_alias_is_supported(tmp_path):
    accepted = write(
        tmp_path / "accepted.json",
        {"status": "source_validated", "release_eligibility": {"source_validation_released": True}},
    )
    rejected = write(
        tmp_path / "rejected.json",
        {"status": "source_validated", "release_flags": {"source_dynamic_validation": "true"}},
    )
    report = filter_results([accepted, rejected], "source_dynamic")
    assert report["accepted_count"] == 1
    assert report["accepted_paths"] == [str(accepted.resolve())]
    assert "release_flag_is_not_boolean" in report["rows"][1]["reasons"]


def test_contradictory_false_gate_overrides_true_alias(tmp_path):
    path = write(
        tmp_path / "conflict.json",
        {
            "release_flags": {"design_eligible": True},
            "gates": {"s015_design_eligible": False},
        },
    )
    row = evaluate(path, "design")
    assert row["accepted"] is False
    assert "contradictory_false_release_flag" in row["reasons"]


def test_numeric_css_result_is_rejected_for_downstream_use(tmp_path):
    path = write(
        tmp_path / "css_only.json",
        {"schema": "result", "status": "css_converged", "product_h2_purity": 0.94},
    )
    assert filter_results([path], "tea_lca")["accepted_count"] == 0


def test_mean_stream_transfer_alias_is_supported_but_still_explicit(tmp_path):
    accepted = write(
        tmp_path / "mean_stream.json",
        {"release_flags": {"dwsim_transfer_allowed": True}},
    )
    rejected = write(
        tmp_path / "closed_mean_stream.json",
        {"release_flags": {"dwsim_transfer_allowed": False}},
    )
    report = filter_results([accepted, rejected], "dwsim")
    assert report["accepted_count"] == 1
    assert report["accepted_paths"] == [str(accepted.resolve())]
