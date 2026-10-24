# PSA-Dyna

**PSA-Dyna: A Physics-Based Dynamic Pressure Swing Adsorption Model**

PSA-Dyna is a Python implementation of a one-dimensional, non-isothermal,
multicomponent pressure swing adsorption (PSA/VPSA) model. It is intended for
transparent research, numerical experiments, independent-case generation and
carefully gated coupling to steady-state process simulators.

> **Research status:** the numerical engine and a figure-level single-bed
> breakthrough benchmark are available. The repository does **not** claim a
> validated industrial cyclic PSA design. Results must pass the explicit
> numerical, source-domain and release gates before they are used in DWSIM,
> process optimization, TEA/LCA or machine-learning datasets.

## Why this repository exists

Steady-state process simulators usually represent PSA as a component splitter.
That approximation cannot reproduce adsorption fronts, pressure equalization,
thermal waves, valve switching or cycle-to-cycle convergence. PSA-Dyna solves
these transient effects in Python and then, only after qualification, converts
a converged cycle to audited cycle-average streams.

The project keeps three questions separate:

1. **Does the numerical model run and conserve mass and energy?**
2. **Does it reproduce an independent experiment within declared limits?**
3. **Is a particular adsorbent, geometry and cycle qualified for process design?**

A positive answer to one question is not treated as evidence for the others.

## Main capabilities

| Capability | Implementation |
|---|---|
| Spatial model | 1D finite-volume packed-bed discretization |
| Gas transport | Advective flux, optional axial dispersion, Ergun pressure-flow closure |
| Adsorption equilibrium | Competitive dual-site Langmuir plus source-local model adapters |
| Adsorption rate | Loading-LDF or pellet partial-pressure-LDF formulations |
| Energy | Gas/solid inventory, adsorption heat, axial conduction, dynamic wall and ambient exchange |
| Gas properties | Ideal gas; optional Peng-Robinson density and departure-energy path |
| Valves | One-way compressible orifice flow with critical choking and near-zero regularization |
| Cycle topology | Independent bed count and step count; four-step and eight-step schedules included |
| Numerics | SciPy stiff ODE integration, sparse Jacobian structure, checkpoint/restart |
| Quality control | Positivity, mass/energy ledgers, full-state CSS, source-domain and release gates |
| Parallel work | Independent-case subprocess execution locally or through an HPC array |
| Process bridge | Guarded cycle-average FEED/PRODUCT/TAIL JSON contract for DWSIM |

See [Physics and numerical formulation](docs/physics.md) for the equations and
implementation choices.

## Installation

Python 3.11 or newer is required.

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e .[dev]
```

Windows PowerShell with an existing Conda environment:

```powershell
& "C:\path\to\python.exe" -m pip install -e ".[dev]"
```

Check the command-line interface:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic --help
```

## Three-minute smoke test

Create a small four-step hydrogen PSA case:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic make-case --source s015 --schedule skarstrom_4step --cells 4 --output "examples\user_s015_4step_4cell.json"
```

Audit the input without solving the ODE system:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic audit-case --case "examples\user_s015_4step_4cell.json" --output "runs\user_s015_input_audit.json"
```

Run one cycle as an execution and bookkeeping test:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic run --case "examples\user_s015_4step_4cell.json" --output "runs\user_s015_smoke" --cycle-budget 1
```

One cycle is a smoke test only. It does not establish cyclic steady state.

## Run to cyclic steady state

Generate a four-bed, eight-step pressure-equalization case:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic make-case --source s015 --schedule equalized_8step --cells 8 --output "examples\user_s015_8step_8cell.json"
```

Run until CSS or the case's `maximum_cycles` limit:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic run --case "examples\user_s015_8step_8cell.json" --output "runs\user_s015_8step_css"
```

Repeat the same command to resume from the last committed checkpoint. A changed
input hash or source-code hash is rejected to prevent mixed-version results.

`--cells 32` means 32 axial finite volumes **per bed**. It does not mean 32
runs, 32 beds or 32 cycle steps.

## Input and output contract

Cases use the versioned `dynamic-psa-case/v1` JSON schema. The important input
groups are:

- species properties;
- bed geometry, voidage, thermal properties and adsorbent layers;
- feed composition, flow, pressure and temperature;
- high, low and product pressure;
- valve areas or source-specific controller data;
- cycle schedule and step durations;
- ODE, trace and CSS tolerances;
- provenance and claim boundary.

Each run produces a self-contained directory with frozen input, materialized
schedule, engineering screen, alternating checkpoints, cycle and step ledgers,
the final cycle trace, a source snapshot and a terminal `result.json` or
structured `failure.json`.

See [Input and output reference](docs/input-output.md).

## Validation status

The current public claim is deliberately narrow.

| Evidence track | Current conclusion |
|---|---|
| Wang et al. H2-rich, single-bed breakthrough | Figure-level comparison passed the declared composition, temperature, mesh, time-step and conservation gates |
| Wang cyclic four-bed VPSA | Not experimentally validated by the single-bed result |
| S-015 four/eight-step process candidate | Numerical design screen only; industrial design gate remains closed |
| Boer et al. CO2/CH4 biogas VPSA | Public figure data were insufficient to identify the apparatus response within the preregistered early-rise gate |
| DWSIM cycle-average transfer | Contract implemented; automatic release currently closed |

For the exact metrics and interpretation, read [Validation evidence and claim
boundaries](docs/validation.md). A machine-readable, copyright-safe summary is
provided in [`validation/validation_summary.json`](validation/validation_summary.json).

## DWSIM interface

The dynamic solver and DWSIM remain separate programs. PSA-Dyna does not copy a
single instantaneous purity into a steady-state flowsheet. It integrates the
converged cycle and prepares three mean streams:

- `FEED`;
- `PRODUCT`;
- `TAIL`.

Each stream carries molar flow, composition, phase-aligned pressure and
temperature, and enthalpy flow. The payload also records the qualified property
package and release blockers. Missing thermodynamic state or evidence closes
the transfer automatically.

See [DWSIM bridge](docs/dwsim-bridge.md).

## Parallel independent cases

PSA-Dyna parallelizes independent cases rather than splitting one small stiff
ODE system across processors:

```powershell
& "C:\path\to\python.exe" -m psa_dynamic.batch --manifest "C:\absolute\path\manifest.json" --workers 2 --require-python 3.11.17
```

Every manifest row has a frozen input SHA-256 and a unique output directory.
One worker cannot overwrite another worker's result.

## Tests

```bash
python -m pytest -q
```

The public suite checks equilibrium derivatives, Peng-Robinson behavior,
vectorized/scalar agreement, schedule topology, input qualification, batch
isolation, restart-related contracts and DWSIM payload rejection behavior.

## Repository layout

```text
PSA-Dyna/
├── src/psa_dynamic/       numerical engine and interfaces
├── examples/              small, cited case definitions
├── tests/                 public regression tests
├── validation/            copyright-safe result summary
├── docs/                  scientific and user documentation
└── .github/workflows/     continuous integration
```

## Documentation

- [Getting started](docs/getting-started.md)
- [Physics and numerical formulation](docs/physics.md)
- [Input and output reference](docs/input-output.md)
- [Validation evidence](docs/validation.md)
- [DWSIM bridge](docs/dwsim-bridge.md)
- [Known limitations](docs/limitations.md)
- [References](docs/references.md)
- [Data and copyright policy](docs/data-policy.md)

## Data and copyright

Publisher PDFs, article figures, proprietary software files and source data
without explicit redistribution permission are intentionally excluded. The
repository cites DOI records and official data repositories and reports only
the derived metrics needed to describe the model's evidence boundary. Users
must obtain any restricted source material from the publisher or author.

## Citation

Please cite the archived software release described in [`CITATION.cff`](CITATION.cff).
The DOI will be added to the main branch after Zenodo completes the first
release archive. No ORCID is asserted by this repository.

## License

The original code and documentation in this repository are released under the
[MIT License](LICENSE). Third-party articles, datasets and trademarks remain
under their respective owners' terms and are not relicensed by this project.
