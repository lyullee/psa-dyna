# Research handoff

This document is the entry point for resuming PSA-Dyna after the original
development session is unavailable. It describes what the repository proves,
what it does not prove, how to reproduce the public state, and which task has
the highest scientific priority.

## One-paragraph project summary

PSA-Dyna is a one-dimensional, non-isothermal, multicomponent dynamic PSA/VPSA
solver written in Python. It was developed to replace an equilibrium splitter
with a transient packed-bed model before coupling hydrogen purification to a
larger anaerobic-digestion-to-hydrogen process study. The software implements
adsorption equilibrium and kinetics, packed-bed mass and energy balances,
pressure-flow closure, valves, cyclic schedules, cyclic-steady-state (CSS)
tests, restartable execution, independent-case parallelism and a guarded
cycle-average DWSIM contract. The public experimental claim is intentionally
narrow: a Wang et al. single-bed H2-rich breakthrough comparison passes its
declared figure-level gates. No public cyclic experiment currently qualifies
the S-015 plant candidate, so design, DWSIM, TEA/LCA, optimization and
training-data release remain closed.

## Authoritative public state

| Item | Value |
|---|---|
| Repository | <https://github.com/lyullee/psa-dyna> |
| Package | <https://pypi.org/project/psa-dyna/> |
| Concept DOI | <https://doi.org/10.5281/zenodo.23608231> |
| v0.1.0 DOI | <https://doi.org/10.5281/zenodo.23608232> |
| Software maturity | Alpha research software |
| Python | 3.11 or newer |
| License | MIT for original repository content |
| Machine-readable status | [`research-state.json`](research-state.json) |
| Validation metrics | [`../validation/validation_summary.json`](../validation/validation_summary.json) |

The Git commit, release tag and archive DOI identify software snapshots. A DOI
does not enlarge the scientific claim boundary.

## What is complete

- A portable NumPy/SciPy dynamic engine exists and is covered by public tests.
- Four-step and eight-step schedules have independent bed-count and step-count
  semantics.
- Run directories freeze inputs, source snapshots, hashes, ledgers,
  checkpoints, traces and terminal status.
- Ordinary-cycle CSS is evaluated from full state changes, not purity alone.
- The public Wang single-bed benchmark passed declared composition,
  temperature, spatial, temporal and conservation screens.
- Independent cases can run in isolated subprocesses locally or through an HPC
  array manifest.
- The mean-stream DWSIM payload and its fail-closed preflight are implemented.
- PyPI packaging, GitHub continuous integration and Zenodo archiving are live.

## What remains scientifically open

1. **Independent cyclic H2 validation.** The single-bed Wang result cannot
   validate the later multibed cycle.
2. **S-015 design qualification.** Numerical CSS and scalar GCI pass, but
   thermal profile, exact material, hydraulic equipment, fixed-target and
   public cyclic gates do not.
3. **Biogas VPSA validation.** The Boer public figures do not identify the IR
   reservoir and apparatus response well enough to pass the preregistered
   early-rise gate.
4. **Independent non-isothermal validation.** A public pressure-temperature
   transient with sufficient geometry, boundary and sensor metadata is still
   needed. CALF-20 and NASA sources are candidates, not released validations.
5. **DWSIM qualification.** The contract exists, but no automatically released
   payload is allowed until the upstream scientific gates pass.

## First hour for a new researcher or agent

Use a fresh directory and do not copy old result folders into it.

```powershell
git clone https://github.com/lyullee/psa-dyna.git
```

```powershell
Set-Location .\psa-dyna
```

```powershell
python -m venv .venv
```

```powershell
& ".\.venv\Scripts\python.exe" -m pip install --upgrade pip
```

```powershell
& ".\.venv\Scripts\python.exe" -m pip install -e ".[dev]"
```

```powershell
& ".\.venv\Scripts\python.exe" -m pytest -q
```

Then read, in order:

1. this handoff;
2. [`research-state.json`](research-state.json);
3. [`decision-log.md`](decision-log.md);
4. [`validation.md`](validation.md);
5. [`known-issues.md`](known-issues.md);
6. [`roadmap.md`](roadmap.md);
7. [`code-map.md`](code-map.md);
8. [`reproducibility-runbook.md`](reproducibility-runbook.md).

## Minimal execution check

Create and audit a small case before spending time on a cyclic run.

```powershell
& ".\.venv\Scripts\python.exe" -m psa_dynamic make-case --source s015 --schedule skarstrom_4step --cells 4 --output "runs\handoff_s015_4cell.json"
```

```powershell
& ".\.venv\Scripts\python.exe" -m psa_dynamic audit-case --case "runs\handoff_s015_4cell.json" --output "runs\handoff_s015_4cell.audit.json"
```

```powershell
& ".\.venv\Scripts\python.exe" -m psa_dynamic run --case "runs\handoff_s015_4cell.json" --output "runs\handoff_s015_smoke" --cycle-budget 1
```

This is an execution check only. It is not a validation or design result.

## Resume rules for long calculations

- Reuse the same case JSON and output directory to resume.
- Never edit the case file or engine source between pauses and resume attempts.
- The runner checks input and engine hashes and rejects mismatches.
- Do not delete `RUNNING.lock` until the owning process is confirmed absent.
- A checkpoint is not a result. Only `result.json` with its stated terminal
  status may be interpreted.
- An accelerated state is only an initialization candidate. Ordinary physical
  cycles must establish the CSS streak after it.
- A paused, failed or unconverged case is never training truth.

## Evidence ladder

The project uses this order. A later level cannot be inferred from an earlier
one.

1. software tests and dimensional checks;
2. positivity and mass/energy conservation;
3. mesh and time-step qualification;
4. ordinary-cycle CSS;
5. fixed-model experimental comparison;
6. exact material and operating-domain qualification;
7. plant thermal, hydraulic and control qualification;
8. downstream DWSIM, TEA/LCA, optimization or training release.

Every downstream consumer should call the release filter or reproduce its
default-deny logic. Numerical convergence alone is never permission.

## Highest-priority next study

Acquire or identify one **independent cyclic H2 PSA dataset** containing, at a
minimum, synchronized step timing, pressure, product/feed flow and composition,
plus bed geometry, adsorbent identity and operating temperature. Freeze its
source locator and hashes, define calibration and no-refit holdouts, and
preregister metrics before running a fit. Do not tune a missing apparatus
volume, sensor lag or valve history against the holdout trace.

If that dataset cannot be obtained, publish the current software and
single-bed validation with an explicitly limited claim. Do not promote S-015
to plant design by adding more numerical resolution alone.

## Public/private boundary

The public repository intentionally excludes publisher PDFs, article figures,
full digitized coordinate sets with uncertain redistribution rights,
proprietary Aspen/DWSIM files, local absolute paths, private process results,
credentials and large checkpoints. A future session should request those
items from the user only when they are essential and should keep them outside
Git unless redistribution permission is documented.

## Updating this handoff

When a scientific gate changes:

1. update `validation/validation_summary.json`;
2. update `docs/research-state.json`;
3. update the relevant validation and issue narrative;
4. add tests for the decision rule;
5. record the change in `CHANGELOG.md`;
6. run the complete test suite;
7. commit the evidence and code together, or record immutable cross-hashes.

