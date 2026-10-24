# Reproducibility runbook

## Environment

The supported minimum is Python 3.11. Use an isolated environment and record
the exact interpreter and dependency versions with every campaign.

```powershell
python --version
```

```powershell
python -m pip install -e ".[dev]"
```

```powershell
python -m pip freeze > environment-lock.txt
```

`environment-lock.txt` is a local run record. Update project dependency bounds
through review rather than committing an arbitrary machine-wide freeze.

## Repository checks

```powershell
git status --short
```

```powershell
git rev-parse HEAD
```

```powershell
python -m pytest -q
```

Do not start a publication campaign from a dirty worktree unless every local
change is intentionally captured in the run's source snapshot.

## Case preparation

Generate a case through the CLI or copy an example to a new file. Never edit a
published example in place during a campaign.

```powershell
python -m psa_dynamic make-case --source s015 --schedule equalized_8step --cells 8 --output "campaign\cases\s015_n8.json"
```

Audit before running:

```powershell
python -m psa_dynamic audit-case --case "campaign\cases\s015_n8.json" --output "campaign\cases\s015_n8.audit.json"
```

Review every `open_engineering_evidence_gates` entry. The audit reports input
evidence; it does not predict convergence.

## One-case execution and resume

```powershell
python -m psa_dynamic run --case "campaign\cases\s015_n8.json" --output "campaign\runs\s015_n8" --cycle-budget 1
```

Run the same command again to continue. Omit `--cycle-budget` for the frozen
case maximum.

Expected run artifacts include:

- `input.json` and `schedule.json`;
- `engineering_input_screen.json`;
- `source_code.zip`;
- alternating `checkpoint.state.*.npz` plus `checkpoint.json`;
- cycle and step ledgers;
- `active_step.json` while solving;
- `result.json` or `failure.json` at a terminal boundary.

Interpret the status before using any scalar value.

## Independent-case campaign

The manifest must contain a frozen engine hash, one SHA-256 per case and a
unique output directory per row. Start with one or two workers.

```powershell
python -m psa_dynamic.batch --manifest "C:\absolute\campaign\manifest.json" --workers 2 --require-python 3.11.17
```

For one array element:

```powershell
python -m psa_dynamic.batch --manifest "C:\absolute\campaign\manifest.json" --index 0 --require-python 3.11.17
```

The dispatcher fixes common BLAS thread counts at one. Each worker is a new
process and writes only its assigned directory.

## Validation workflow

1. Freeze source DOI, exact locator, license/access note and file hashes.
2. Freeze units and conversions.
3. Separate calibration and holdout observations.
4. Freeze metrics, alignment rule and thresholds.
5. Run software and conservation checks.
6. Run spatial and temporal refinement.
7. Establish ordinary CSS for cyclic cases.
8. Fit only declared calibration parameters.
9. Run holdouts once without refitting.
10. Write explicit boolean release flags and the claim boundary.

See [`validation-playbook.md`](validation-playbook.md) for details.

## Downstream release check

The filter accepts only an explicit boolean true for the requested purpose.
For example, a future qualified result can be checked with:

```powershell
python -m psa_dynamic.release_filter --purpose dwsim --input "campaign\runs\qualified\result.json" --output "campaign\release\dwsim_gate.json" --require-nonempty
```

A non-zero exit with zero accepted records is a correct fail-closed result.

## Publication archive checklist

- clean commit and immutable tag;
- passing public tests;
- case, code and result hashes;
- interpreter and dependency versions;
- calibration/holdout declaration;
- mesh, time-step, CSS and balance evidence;
- exact claim boundary;
- copyright review;
- updated `validation_summary.json` and `research-state.json`;
- archive DOI linked to the exact release.

