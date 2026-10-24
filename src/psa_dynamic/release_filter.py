"""Fail-closed selection of PSA results for downstream use.

The dynamic engine produces many useful diagnostics.  A completed run, a
small residual, or a passing numerical screen is *not* permission to use that
run as design, DWSIM, TEA/LCA, source-validation, or training data.  This
module is the single small gate that downstream scripts can call before
reading a result file.

Only an explicit boolean ``true`` release flag for the requested purpose is
accepted.  Missing, false, string-valued, or conflicting flags are rejected.
The default is therefore safe when a new result schema is introduced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


PURPOSES = ("source_dynamic", "design", "dwsim", "tea_lca", "training")

# Several existing audits use slightly different names.  Supporting those
# aliases here is intentional, but every alias is still required to be the
# literal JSON boolean true.  No truthy strings or numeric values are allowed.
FLAG_ALIASES: dict[str, tuple[str, ...]] = {
    "source_dynamic": ("source_dynamic_validation", "source_validation_released"),
    "design": ("design_eligible", "s015_design_eligible"),
    "dwsim": ("dwsim_application_allowed", "dwsim_transfer_allowed", "dwsim_mean_stream"),
    "tea_lca": ("tea_lca_allowed", "tea_lca_transfer_allowed", "tea_lca"),
    "training": ("training_eligible",),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _flag_maps(document: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    """Return all supported flag maps, preserving their provenance."""

    maps: list[tuple[str, dict[str, Any]]] = []
    for key in ("release_flags", "release_eligibility", "gates"):
        value = document.get(key)
        if isinstance(value, dict):
            maps.append((key, value))
    # A few compact result records expose flags at the top level.  Include the
    # top-level record last, while avoiding a duplicate nested map.
    maps.append(("top_level", document))
    return maps


def _resolve_flag(document: dict[str, Any], purpose: str) -> tuple[str | None, Any]:
    aliases = FLAG_ALIASES[purpose]
    for map_name, mapping in _flag_maps(document):
        for alias in aliases:
            if alias in mapping:
                return f"{map_name}.{alias}", mapping[alias]
    return None, None


def evaluate(path: Path, purpose: str) -> dict[str, Any]:
    """Evaluate one JSON result and return an auditable decision."""

    path = Path(path).resolve()
    row: dict[str, Any] = {
        "path": str(path),
        "purpose": purpose,
        "exists": path.exists(),
        "sha256": None,
        "accepted": False,
        "reasons": [],
    }
    if not path.exists():
        row["reasons"].append("file_missing")
        return row
    row["sha256"] = sha256(path)
    try:
        document = json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception as exc:  # pragma: no cover - defensive CLI boundary
        row["reasons"].append(f"invalid_json:{type(exc).__name__}")
        return row
    if not isinstance(document, dict):
        row["reasons"].append("root_is_not_object")
        return row
    row["status"] = document.get("status")
    row["schema"] = document.get("schema")
    location, value = _resolve_flag(document, purpose)
    row["flag_location"] = location
    row["flag_value"] = value
    if location is None:
        row["reasons"].append("required_release_flag_missing")
    elif type(value) is not bool:
        row["reasons"].append("release_flag_is_not_boolean")
    elif value is not True:
        row["reasons"].append("release_flag_false")
    else:
        row["accepted"] = True
    # A contradictory explicit false flag always wins.  This protects against
    # a document containing both an alias in one map and a false gate in
    # another map.
    false_aliases = set(FLAG_ALIASES[purpose])
    explicit_false = []
    for map_name, mapping in _flag_maps(document):
        for alias in false_aliases:
            if alias in mapping and type(mapping[alias]) is bool and mapping[alias] is False:
                explicit_false.append(f"{map_name}.{alias}")
    if explicit_false:
        row["contradictory_false_flags"] = explicit_false
        row["accepted"] = False
        row["reasons"].append("contradictory_false_release_flag")
    return row


def filter_results(paths: Iterable[Path], purpose: str) -> dict[str, Any]:
    if purpose not in PURPOSES:
        raise ValueError(f"unknown purpose: {purpose}; choose one of {PURPOSES}")
    rows = [evaluate(Path(path), purpose) for path in paths]
    accepted = [row["path"] for row in rows if row["accepted"]]
    rejected = [row["path"] for row in rows if not row["accepted"]]
    return {
        "schema": "dynamic-psa-release-filter/v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": purpose,
        "policy": "default_deny_explicit_boolean_release_flag_required",
        "input_count": len(rows),
        "accepted_count": len(accepted),
        "rejected_count": len(rejected),
        "accepted_paths": accepted,
        "rejected_paths": rejected,
        "rows": rows,
        "interpretation": (
            "Only results with an explicit JSON boolean true release flag for "
            "the requested purpose are returned. Numerical convergence alone "
            "cannot open a downstream release gate."
        ),
    }


def _expand_inputs(inputs: list[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        if item.is_dir():
            paths.extend(sorted(item.rglob("*.json")))
        else:
            paths.append(item)
    # Keep the report deterministic and avoid duplicate paths from overlapping
    # directory/file arguments.
    return list(dict.fromkeys(path.resolve() for path in paths))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--purpose", choices=PURPOSES, required=True)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--require-nonempty",
        action="store_true",
        help="return a non-zero exit code when no result passes the gate",
    )
    args = parser.parse_args()
    report = filter_results(_expand_inputs(args.input), args.purpose)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("schema", "purpose", "input_count", "accepted_count", "rejected_count", "output") if k != "output"}, ensure_ascii=False))
    print(json.dumps({"output": str(args.output)}, ensure_ascii=False))
    return 2 if args.require_nonempty and report["accepted_count"] == 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
