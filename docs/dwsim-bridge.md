# DWSIM cycle-average bridge

## Purpose

DWSIM is a steady-state process simulator. PSA-Dyna is a transient cyclic
solver. The bridge converts a qualified CSS cycle into audited average streams
without pretending that an instantaneous state is a steady-state product.

## Payload schema

The nested payload uses `dynamic-psa-dwsim-mean-stream-input/v1` and contains
the common species order plus `FEED`, `PRODUCT` and `TAIL` streams. Each stream
contains:

- `molar_flow_mol_s`;
- `mole_fractions`;
- `pressure_bar`;
- `temperature_k`;
- `enthalpy_flow_w`;
- a description of the state basis.

The payload also records units, numerical qualification, the property-package
reference, blockers and release flags.

## Build a mean-stream report

```text
python -m psa_dynamic.mean_stream_bridge \
  --result RUN/result.json \
  --case CASE.json \
  --output BRIDGE.json \
  --spatial-audit SPATIAL.json \
  --temporal-audit TEMPORAL.json \
  --profile-audit PROFILE.json
```

All three numerical audits must be supplied together. The bridge also requires
ordinary CSS, the required CSS streak and balance residuals within the frozen
case tolerances.

## Structural preflight

```text
python -m psa_dynamic.dwsim_payload_preflight BRIDGE.json --output PREFLIGHT.json
```

The preflight accepts either the complete bridge report or the nested DWSIM
payload. It checks component order, positive flow, normalized composition,
finite enthalpy, resolved pressure/temperature and a complete property-package
reference.

## Fail-closed policy

A cycle-average ledger alone does not define one scalar pressure or
temperature. PSA-Dyna therefore does not substitute the high-pressure boundary
as a fictitious mean product state. Phase-aligned stream state is required.

The property-package record must identify its name, version, source and SHA-256.
An optional file path can be checked against the hash. Structural preflight
never authorizes transfer by itself.

At version 0.1.0, the JSON bridge and preflight are implemented. Direct
mutation of a DWSIM flowsheet is intentionally outside the public release.
