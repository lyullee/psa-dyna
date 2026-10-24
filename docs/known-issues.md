# Known issues and failed approaches

This page records scientific and operational problems that a future session
must not rediscover from scratch.

## Scientific blockers

### K-001: single-bed validation is not cyclic validation

**State:** open.

The Wang H2-rich breakthrough comparison passed, but it does not exercise
multibed synchronization, blowdown, purge, equalization, product receiver
dynamics or cyclic recovery. Do not use that pass to label a cyclic design as
validated.

**Resolution needed:** independent synchronized cyclic data and a no-refit
holdout comparison.

### K-002: S-015 design release is closed

**State:** open.

N512 ordinary CSS and scalar GCI pass. The 256-to-512 phase-aligned gas and
wall temperature RMS differences were approximately 1.3468 K and 1.4829 K,
above the frozen 1 K limits. Exact material grade/binder, plant heat transfer,
compressor/pipe/valve behavior and public cyclic validation are also missing.
The product result was 101.9625 kg/h H2 versus a 103.8761 kg/h target, with
0.889756 purity.

**Do not:** widen limits, fit equipment after seeing the target, or interpret
more cells as material qualification.

### K-003: Boer early-rise response is not identifiable

**State:** blocked by public source information.

The 18-case apparatus screen achieved a best normalized-profile RMSE of about
0.00832, but the best relative `t05` error was about 8.58%, above the frozen 5%
limit. The source identifies an IR-reservoir contribution but not the reservoir
volume or raw response history.

**Do not:** introduce an arbitrary reservoir volume, detector lag or time shift
after inspecting the validation trace.

**Resolution needed:** raw apparatus response or independently measured
reservoir and sensor dynamics.

### K-004: non-isothermal independent validation is incomplete

**State:** open.

CALF-20 and NASA sources contain useful candidates, but public availability of
a figure, workbook or simulation does not automatically provide synchronized
pressure, temperature, boundary and sensor metadata.

**Resolution needed:** a source-complete independent transient and a dedicated
adapter whose parameters never leak into other branches.

### K-005: wet and contaminated gas are outside the dry PSA claim

**State:** open.

Water, H2S and trace contaminants can dominate guard-bed life and regeneration.
The present public examples must not be treated as wet-biogas equipment
designs.

## Numerical and operational issues

### K-006: CSS can look converged in scalar outputs while state drifts

**Symptom:** purity and recovery change little, but trace composition, loading
or temperature controls the residual.

**Practice:** use `state_scaled_max`, its location and a consecutive ordinary
cycle streak. Inspect both scalar outputs and phase-aligned profiles.

### K-007: accelerated states can be falsely reassuring

**Symptom:** an extrapolated state has a smaller apparent residual, followed by
growth on the next physical map.

**Practice:** retain the proposal as a candidate only, apply growth and
conservation guards, then require ordinary-cycle postchecks. Rejected proposals
must not contaminate the accepted history.

### K-008: long stiff runs may appear frozen

**Symptom:** little console output during a BDF step.

**Practice:** inspect `active_step.json`, process ownership and the worker log.
Respect `maximum_step_wall_s`. A structured wall-budget failure is preferable
to deleting locks or silently restarting.

### K-009: resume can mix incompatible states if hashes are bypassed

**Symptom:** a checkpoint exists after code or input changes.

**Practice:** never bypass input/code/state hash checks. Start a new output
directory for a changed case. Preserve the old terminal record.

### K-010: pressure and product-port semantics can be misread

**Symptom:** a closed or receiving bed-end cell is reported as a product
outlet, or partial pressure is interpreted as total pressure.

**Practice:** use the materialized schedule and phase-aligned boundary ledger.
Score the actual product path, not a convenient terminal cell.

### K-011: more local workers can make runs slower

**Symptom:** CPU oversubscription and memory pressure.

**Practice:** parallelize independent cases, set BLAS thread counts to one,
and increase workers only after measuring host memory and wall time.

## Integration issues

### K-012: a DWSIM payload can be numerically complete but scientifically closed

**State:** expected behavior.

The bridge checks thermodynamic fields and evidence flags. Missing or false
release flags reject transfer even if stream numbers exist.

### K-013: property packages are not interchangeable labels

Peng-Robinson in Python/DWSIM and Aspen PR-MHV2 do not become equivalent by
using similar names. A property-method comparison must document components,
binary parameters, state range and reference enthalpy.

## Documentation and provenance issues

### K-014: public figures are not raw historians

Digitized curves have axis, raster, timing and synchronization uncertainty.
Call them figure-level comparisons and retain the digitization protocol. Do
not claim raw-data equivalence.

### K-015: copyrighted and proprietary artifacts cannot be committed casually

Keep publisher PDFs, figures, supplementary files without redistribution
permission, Aspen backups, DWSIM binaries and private integrated-process files
outside the public repository. Cite and hash them locally.

## Reporting a new issue

Record:

- exact commit and Python/SciPy versions;
- case and engine hashes;
- command line;
- terminal status and failure JSON;
- whether the result was calibration, holdout, numerical or design evidence;
- which release flags remained false;
- the smallest reproducible case that preserves the failure.

