# Development history

This is a public-safe reconstruction of the development path. It records the
scientific decisions and reversals needed to understand the current code. It
is not a verbatim chat log, and it does not expose private paths or restricted
source files.

## Phase 1: process-model motivation

The wider research project connects anaerobic digestion, biogas conditioning,
reforming, water-gas shift, hydrogen purification and TEA/LCA. DWSIM is the
primary open process simulator and Aspen Plus is used for cross-checking.

The first purification representation was necessarily simple: a separator
with prescribed split fractions. That is useful for flowsheet bookkeeping but
cannot answer questions about bed count, step timing, pressure equalization,
temperature fronts, cyclic convergence or robustness. The PSA work therefore
became a separate physics-model project.

## Phase 2: common dynamic engine

The model was built as a one-dimensional finite-volume packed bed with:

- multicomponent gas and adsorbed-phase inventories;
- competitive adsorption equilibrium;
- LDF mass transfer;
- gas/solid and wall energy storage;
- heats of adsorption and heat exchange;
- Ergun pressure-flow closure;
- compressible one-way valves;
- configurable multi-bed cycle schedules.

Bed count and cycle-step count were deliberately separated. A four-step cycle
does not mean four beds, and an eight-step schedule is not represented by
running a four-step model twice.

## Phase 3: numerical hardening

Dynamic PSA is stiff and cyclic. Several practices were introduced after early
runs proved that a final purity value was not enough:

- full-state CSS residuals replaced scalar-only convergence claims;
- CSS requires consecutive ordinary physical cycles;
- mass and energy ledgers are integrated over each step and cycle;
- alternating atomic checkpoints protect against interrupted writes;
- input and engine hashes prevent mixed-version resume;
- a source-code snapshot is stored with each run;
- wall-time guards create structured failures instead of silent hangs;
- acceleration proposals must be followed by ordinary-cycle postchecks;
- spatial and temporal refinement are separate gates.

Independent cases are dispatched as separate processes. This is more robust
than trying to parallelize one small stiff ODE solve and maps naturally to HPC
job arrays.

## Phase 4: hydrogen evidence

Wang et al. was split into two distinct questions:

1. can the engine reproduce a single packed-bed H2-rich breakthrough and its
   thermal response?;
2. can it reproduce and design a cyclic multibed VPSA?

The first question passed the declared public-figure, mesh, time-step and
conservation gates without parameter fitting in the validation workflow. The
second question remains open because a public synchronized cyclic historian is
not available in the present evidence package.

This separation prevented a common error: treating a breakthrough fit as proof
of cyclic purity, recovery or plant capacity.

## Phase 5: S-015 engineering candidate

An S-015 layered activated-carbon/5A process candidate was created for the
larger hydrogen flowsheet. Long runs and safeguarded numerical initialization
eventually established an ordinary N512 periodic state and a passing scalar
128/256/512 GCI screen.

That numerical success did not close the design gate. The phase-aligned
thermal profiles exceeded the frozen 1 K comparison limits, an ordinary N1024
CSS result was not established, exact commercial adsorbent equivalence was not
shown, and plant heat-transfer and equipment behavior were not qualified. The
fixed product calculation also returned 101.9625 kg/h H2 against a 103.8761
kg/h target and 0.889756 H2 purity.

The lesson is explicit: increasing mesh resolution cannot supply missing
material, equipment or experimental evidence.

## Phase 6: biogas and CCS evidence search

The validation search expanded beyond hydrogen to CO2/CH4 biogas VPSA and
non-isothermal CCS datasets. Each source received its own adapter and parameter
namespace. Parameters were never pooled across gases or adsorbents.

The Boer branch was especially instructive. The overall normalized profile
could look close while the early-rise timing failed. Eighteen preregistered
apparatus-decomposition cases produced a best profile RMSE near 0.00832, but
the minimum relative `t05` error was about 8.58%, above the 5% limit. Because
the public source did not identify the IR-reservoir volume or raw apparatus
response, the missing response was not tuned after the fact. The branch was
closed as a negative identifiability result.

CALF-20 and NASA material were retained as independent non-isothermal
candidates. They were not re-labelled as validation merely because figures or
computational repositories existed.

## Phase 7: downstream release control

A cycle-average bridge was implemented for DWSIM. It transfers mean FEED,
PRODUCT and TAIL states only after checking phase-aligned flow, composition,
pressure, temperature, enthalpy reference and evidence flags.

The `release_filter` module introduced a default-deny policy. Missing flags,
false flags, string-valued flags and contradictory flags all reject a result.
This prevents a successful numerical run from silently entering DWSIM,
TEA/LCA, optimization or machine-learning data.

## Phase 8: public release

The reusable engine, small examples, tests, derived validation metrics and
documentation were separated from the private integrated-process workspace.
Publisher PDFs, article figures, proprietary simulator files, large run
folders and private paths were excluded.

Version 0.1.0 was published through GitHub, PyPI and Zenodo. The release is an
auditable software snapshot, not an industrial-design certificate.

## Current position

The code-development objective for a reusable physics engine is substantially
complete. The main remaining work is evidence acquisition and independent
validation. The most valuable next result is not another optimizer or finer
mesh; it is a source-complete cyclic H2 dataset and a preregistered no-refit
comparison.

