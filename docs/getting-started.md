# Getting started

## Requirements

- 64-bit Python 3.11 or newer
- NumPy
- SciPy
- pytest for development and regression tests

DWSIM, Aspen Plus, COM, MPI and a GPU are not needed to run the Python model.

## Recommended workflow

1. Generate or copy a case JSON.
2. Run `audit-case` before an expensive simulation.
3. Perform a one-cycle smoke test on a coarse mesh.
4. Inspect `failure.json`, balances and positivity checks.
5. Run to ordinary CSS.
6. Repeat on a refined spatial grid and smaller maximum time step.
7. Apply the source-domain and release audits before downstream use.

## Command reference

### Generate a case

```text
python -m psa_dynamic make-case \
  --source {wang,s015,guderian} \
  --schedule SCHEDULE \
  --cells N \
  --output CASE.json
```

The Guderian adapter accepts experiment IDs 13, 14 and 15 with
`--experiment-id` and requires the `guderian_two_bed_6step_source` schedule.

### Audit a case

```text
python -m psa_dynamic audit-case --case CASE.json --output AUDIT.json
```

This is a read-only engineering screen. It identifies missing evidence and
out-of-domain parameters without solving the dynamic system.

### Run or resume

```text
python -m psa_dynamic run --case CASE.json --output RUN_DIRECTORY
```

Optional controls:

- `--cycle-budget N`: finish at most N additional complete cycles;
- `--step-budget N`: finish at most N additional physical schedule slots.

The output directory is single-writer. Never run two processes against the
same directory.

## Reading a completed run

Start with `result.json`. Confirm all of the following before interpreting
purity or recovery:

1. `status` is the expected terminal state;
2. `state_css_converged` is true for a cyclic case;
3. mass and energy residuals are within the frozen case limits;
4. `trajectory_domain_screen.json` reports no material-domain violation;
5. the requested downstream release flag is true.

A numerically completed result can still be unsuitable for design.
