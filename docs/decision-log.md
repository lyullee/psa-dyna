# Decision log

The decisions below explain why the repository has its present structure.
They should be changed only with a documented scientific reason.

## D-001: physics model before surrogate model

**Decision:** build and qualify a conservation-based dynamic model before
generating data for a physics-informed or other machine-learning surrogate.

**Reason:** a surrogate can reproduce systematic model error. Unqualified,
paused or unconverged cases would create a large but misleading dataset.

**Consequence:** `training_eligible` defaults to false.

## D-002: one-dimensional finite volumes

**Decision:** use a 1D axial finite-volume packed-bed model.

**Reason:** it captures the dominant axial fronts and cyclic interactions at a
cost suitable for case campaigns.

**Boundary:** it cannot resolve radial gradients, channeling, maldistribution,
distributor details, attrition or pellet cracking.

## D-003: source-isolated parameter namespaces

**Decision:** keep every paper, gas system and adsorbent in a separate adapter.

**Reason:** equal chemical names do not imply equal adsorbent grade, binder,
pore structure, fitted domain, kinetics or heat data.

**Consequence:** no Wang, Ahn, Boer, CALF-20, NASA or other parameters may be
mixed unless a new, explicit cross-source model and validation justify it.

## D-004: bed count and step count are independent

**Decision:** schedule topology represents beds, steps and connections
explicitly.

**Reason:** “four-step” and “four-bed” answer different questions. Pressure
equalization requires donor/receiver semantics that a repeated separator
cannot reproduce.

## D-005: full-state CSS

**Decision:** CSS is based on scaled changes in the full process state across
ordinary cycles and requires a consecutive passing streak.

**Reason:** purity can appear stationary while internal loading, trace species
or temperature continues to drift.

## D-006: acceleration is not evidence

**Decision:** an accelerated image is an initial guess only. Ordinary physical
maps must pass the unchanged CSS and conservation limits afterwards.

**Reason:** an accelerator can reduce an algebraic residual while leaving the
physical periodic trajectory or a trace component unconverged.

## D-007: separate numerical and experimental gates

**Decision:** conservation, mesh refinement, time refinement, CSS,
experimental comparison and design qualification are separate gates.

**Reason:** each detects a different failure. Passing one cannot logically
prove another.

## D-008: preregister calibration and holdout use

**Decision:** declare fitted parameters, calibration traces, holdouts, metrics
and thresholds before examining holdout results.

**Reason:** post-hoc boundary, sensor or apparatus fitting can turn missing
metadata into an apparently excellent but non-identifiable match.

## D-009: preserve negative results

**Decision:** retain failed gates and structured failures.

**Reason:** the Boer early-rise miss and other terminal failures define what
the current model and public evidence cannot support. Hiding them would cause
future sessions to repeat work or overstate validation.

## D-010: ideal gas by default, Peng-Robinson as an option

**Decision:** make the gas-property closure explicit in each case.

**Reason:** most validation ranges can be screened with an ideal model, while
denser states may need Peng-Robinson. Selecting PR does not eliminate the need
for component and binary-interaction qualification.

## D-011: case-level parallelism

**Decision:** parallelize independent immutable cases as subprocesses or HPC
array jobs; force one BLAS thread per worker.

**Reason:** it gives clean failure isolation and usually scales better than
splitting one stiff, modest-size ODE system.

## D-012: hash-bound restart

**Decision:** resume only when case, engine and checkpoint hashes match.

**Reason:** silently continuing an old state with changed equations or input
creates a trajectory that belongs to no reproducible model.

## D-013: default-deny downstream release

**Decision:** require an explicit JSON boolean true for each downstream
purpose.

**Reason:** successful execution and even numerical convergence do not qualify
a result for design, DWSIM, TEA/LCA, optimization or training.

## D-014: cycle-average DWSIM coupling

**Decision:** transfer integrated mean FEED, PRODUCT and TAIL states rather
than a single instantaneous purity.

**Reason:** DWSIM is steady state; the only defensible reduced interface is an
audited cycle-average material and energy contract.

## D-015: copyright-safe public repository

**Decision:** publish original code, tests, small cited inputs and derived
metrics, but not publisher figures, restricted source data or proprietary
simulator files.

**Reason:** reproducibility must respect redistribution rights. Source hashes
and precise locators can establish provenance without re-publishing protected
content.

