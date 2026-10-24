# Validation playbook

## Purpose

This playbook prevents three common substitutions:

- code verification for experimental validation;
- a breakthrough curve for cyclic performance;
- a converged simulation for plant qualification.

## Evidence levels

| Level | Question | Minimum evidence |
|---|---|---|
| L0 | Can the case be parsed? | Schema, units and topology checks |
| L1 | Is the implementation internally consistent? | Tests, positivity and ledgers |
| L2 | Is the numerical answer resolved? | Mesh, time-step and CSS evidence |
| L3 | Does the fixed model reproduce an experiment? | Preregistered, no-refit holdouts |
| L4 | Is the material domain applicable? | Exact grade/binder and mixture domain |
| L5 | Is plant use defensible? | Thermal, hydraulic, valve/control and safety evidence |
| L6 | May a downstream tool consume it? | Explicit purpose-specific release flag |

## Source contract

For every validation adapter, record:

- DOI and exact table, figure or repository locator;
- file SHA-256 where a legally obtained local artifact is used;
- redistribution license or the reason the artifact stays private;
- gas species and basis;
- adsorbent name, grade, binder and batch when known;
- geometry, voidage and packing data;
- temperature, total and partial-pressure domain;
- boundary histories and sensor locations;
- all unit conversions;
- which parameters are published, inferred, fixed or fitted.

An unresolved item is a gate, not an invitation to borrow a value from another
paper.

## Calibration and holdout protocol

Before executing the validation:

1. name the calibration observations;
2. name every holdout observation;
3. list adjustable parameters and bounds;
4. state the alignment and interpolation rule;
5. set each acceptance limit;
6. freeze code and input hashes.

After a holdout is viewed, changing these items creates a new exploratory
analysis. It cannot replace the original confirmatory result silently.

## Required numerical checks

### Conservation

Report relative component and total mass closure and energy closure over the
same interval used for performance metrics. A small solver residual is not a
balance.

### Spatial refinement

Compare both scalar cycle metrics and phase-aligned profiles. A passing purity
change can coexist with an unresolved thermal front.

### Temporal refinement

Reduce the maximum ODE step or otherwise tighten the temporal control while
preserving the model and input. Report the differences against frozen limits.

### CSS

Use the full scaled process state across ordinary cycles. Require consecutive
passing cycles and report the controlling variable, bed and cell. Warm starts
and accelerated images do not inherit CSS.

## Experimental scoring

Prefer synchronized pressure, temperature, composition and flow. Report
profile and landmark metrics because one can hide the other's failure.

Useful metrics include:

- composition RMSE on a declared grid;
- temperature RMSE and peak location/time;
- `t05`, `t50` and `t95` breakthrough landmarks;
- pressure-step rise and decay timing;
- cycle-integrated feed, product and tail amounts;
- purity, recovery and productivity;
- repeatability or measurement uncertainty.

Do not time-shift, filter or add sensor lag unless the transformation was
declared or independently measured.

## Parameter isolation

The same numerical engine may host multiple gases and adsorbents. That does
not authorize common parameters. Each adapter must materialize its own source
contract. Cross-source transfer is a separate hypothesis requiring its own
test.

## Decision record

A validation report should end with machine-readable booleans for:

- source dynamic validation;
- design eligibility;
- DWSIM transfer;
- TEA/LCA transfer;
- optimization;
- training-data use.

All default to false. State the narrowest supported claim in plain language.

## Current examples

- **Wang single-bed:** L3 figure-level pass for the declared breakthrough
  observables; no cyclic or design release.
- **S-015:** selected L2 numerical passes; L3-L5 incomplete or failed.
- **Boer:** negative L3 apparatus-identifiability result under the frozen
  early-rise gate.

