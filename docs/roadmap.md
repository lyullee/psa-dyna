# Research roadmap

The order below maximizes scientific value. It is intentionally evidence-first
rather than feature-first.

## Priority 1: independent cyclic hydrogen validation

**Goal:** validate a complete cyclic H2 separation case, not only a
breakthrough column.

**Required source fields:** synchronized step events, bed pressure,
feed/product flow, product composition, geometry, adsorbent identity,
temperature and enough valve/control information to reproduce boundaries.

**Deliverables:** source contract, frozen calibration/holdout split, adapter,
mesh/time/CSS package, no-refit holdout report and explicit release flags.

**Stop condition:** if a decisive boundary or apparatus response is absent,
publish the gap rather than fitting it to the holdout.

## Priority 2: S-015 qualification

Proceed only after Priority 1 supplies a credible cyclic basis.

- close exact activated-carbon and 5A grade/binder evidence;
- qualify competitive mixture range over the accepted trajectory;
- resolve phase-aligned thermal profile convergence;
- establish plant radial heat-transfer evidence;
- qualify compressor, pipe and control-valve behavior;
- meet the fixed H2 product, purity and recovery envelope;
- repeat the release audit without weakening frozen limits.

## Priority 3: CO2/CH4 biogas VPSA

Seek raw or openly licensed apparatus-response data for the Boer system, or a
different source-complete biogas VPSA dataset. Preserve the current 18-case
negative result. Do not overwrite it with a new sensor/reservoir fit.

## Priority 4: independent non-isothermal transient

Evaluate CALF-20 and NASA candidates under separate namespaces. Select a case
only if pressure, temperature, geometry, boundary and sensor semantics are
sufficient. A computational repository may be a code cross-check but is not an
experimental holdout by itself.

## Priority 5: DWSIM and integrated process use

After the upstream gates pass:

1. validate the thermodynamic state and enthalpy reference;
2. produce cycle-average FEED/PRODUCT/TAIL streams;
3. run payload preflight;
4. compare the DWSIM steady-state material and energy balances;
5. then permit TEA/LCA and flowsheet optimization.

No current result is automatically released for this step.

## Priority 6: surrogate and large campaigns

Only qualified cases may become labels. Start with a bounded design of
experiments, keep train/validation/test regions separate, preserve failure
classes and run independent cases through local subprocesses before migrating
the same manifest to an HPC array.

## Publication options

### Conservative software/validation paper

Publish the engine, numerical qualification, Wang single-bed figure-level
benchmark, negative Boer identifiability result and fail-closed coupling
framework. The claim must exclude industrial cyclic validation.

### Stronger process-design paper

Requires a passed independent cyclic H2 validation and closed S-015 material,
thermal, hydraulic and target-performance gates before DWSIM/TEA conclusions.

## Work not worth prioritizing now

- another optimization algorithm on unqualified design outputs;
- larger mesh runs that do not address missing source/equipment evidence;
- AI training on unconverged or unreleased cases;
- cosmetic agreement obtained by post-hoc time shifts or apparatus fitting.

