# Validation evidence and claim boundaries

## Validation philosophy

PSA-Dyna separates software verification, numerical qualification,
experimental comparison and plant-design qualification.

- **Software verification** asks whether code paths and ledgers behave as
  intended.
- **Numerical qualification** asks whether mesh and time-step refinement are
  adequate and CSS is genuine.
- **Experimental comparison** asks whether a fixed model reproduces measured
  observables within preregistered limits.
- **Design qualification** additionally requires the correct commercial
  adsorbent, equipment, thermal, valve and control information for the target
  plant.

## Wang 2023 H2-rich single-bed breakthrough

Source: Wang et al., *Chemical Engineering Journal* 472 (2023) 144920,
DOI `10.1016/j.cej.2023.144920`.

The public comparison uses the article's single packed-bed breakthrough
experiment with activated carbon RB3:

- pressure: 5 bar;
- temperature: 293.95 K;
- feed: 58% H2, 29% CH4, 7% N2 and 6% CO2;
- observables: outlet H2/CH4/N2/CO2 and five axial temperatures.

No parameter was fitted in the current validation workflow. Meshes of 16, 32
and 64 axial cells were run at an 8 s maximum ODE step; the 64-cell case was
repeated at 4 s.

### Numerical results

| Test | Observed value | Declared limit | Result |
|---|---:|---:|---|
| 32-to-64 cell maximum outlet RMSE | 0.003523 | 0.01 | pass |
| 32-to-64 cell maximum temperature RMSE | 0.302646 K | 0.5 K | pass |
| 8-to-4 s maximum outlet RMSE | 1.44e-8 | 0.01 | pass |
| 8-to-4 s maximum temperature RMSE | 1.92e-6 K | 0.5 K | pass |
| Worst reported component-balance residual | < 2e-9 | 1e-5 | pass |
| Worst reported energy-balance residual | < 1.6e-7 | 1e-5 | pass |

The declared figure-level composition and temperature gates passed.

### What this result does not prove

The comparison used digitized public-figure information, not a raw synchronized
historian. It validates neither the four-bed cyclic VPSA experiment nor the
S-015 process design. Raw-source equivalence, cyclic validation, DWSIM release,
TEA/LCA release and training-data release remain false.

## S-015 process candidate

The S-015 case is a literature-parameterized layered activated-carbon/5A
engineering screen. It is not the Wang apparatus.

The current high-resolution branch reached ordinary CSS at 512 cells and
passed a scalar 128/256/512 GCI screen. Nevertheless, the design gate remains
closed:

- predicted H2 product: 101.9625 kg/h;
- fixed target: 103.8761 kg/h;
- target error: -1.842%;
- predicted H2 purity: 0.889756.

Open evidence includes exact adsorbent grade/binder and mixture fit domain,
plant-scale heat transfer, selected compressor/pipe/control-valve behavior,
fixed-target product performance and independent cyclic experimental data.
Additional mesh refinement alone cannot close those gaps.

A minimum literature-domain supplement now maps the accepted N512 5A-layer
trajectory (293.15--368.59 K and 1.10--15.14 bar) against source-specific
experimental scopes. The closest direct quaternary equilibrium source covers
298.15--318.15 K and 1--18 bar but does not establish Ahn-grade equivalence.
The 293.15--298.15 K fringe has only pure-component support, and the hotter
part of the trajectory remains an explicit extrapolation.
The implemented bed-to-wall value is exactly `k_eff/R`, not an independently
qualified wall coefficient, and the frozen N256-to-N512 gas and wall
temperature RMS limits are not met. See [S-015 minimum applicability
supplement](s015-minimum-applicability.md).

## CO2/CH4 biogas VPSA evidence

The Boer public-figure branch was handled as a negative identifiability result.
Eighteen preregistered apparatus-decomposition cases were completed. The best
profile RMSE was approximately 0.00832 in normalized concentration, but the
minimum early-rise (`t05`) relative error was approximately 8.58%, above the
predeclared 5% limit.

The source attributes an important tail to the IR reservoir but does not
publish the reservoir volume or a raw synchronized response. The project did
not add a post-hoc detector or reservoir parameter. This branch therefore
remains closed pending raw response data or independently measured apparatus
volumes and sensor dynamics.

## Automated release behavior

A run may be numerically successful and still remain closed for downstream
use. Release filters inspect source validation, material-domain coverage,
spatial/time/profile qualification, phase-aligned thermodynamic state and the
named DWSIM property-package reference.

The machine-readable summary in `validation/validation_summary.json` contains
derived metrics only. It deliberately excludes publisher figures, PDFs and
restricted raw data.
