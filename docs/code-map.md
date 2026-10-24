# Code and artifact map

This map helps a new developer locate the scientific responsibility of each
public file without reconstructing the architecture from imports.

## Execution path

```text
case JSON
  -> cases.py / source adapter
  -> schedules.py
  -> qualification.py input screen
  -> core.py Network and packed-bed equations
  -> numerics.py Jacobian and quadrature helpers
  -> runner.py BDF integration, ledgers, CSS and checkpointing
  -> result.json or failure.json
  -> release_filter.py
  -> mean_stream_bridge.py
  -> dwsim_payload_preflight.py
```

The release filter is intentionally downstream of the solver. The solver may
produce scientifically useful negative or diagnostic results without granting
permission to use them in another model.

## Package modules

| Module | Responsibility | Change risk |
|---|---|---|
| `__main__.py` | User CLI for case creation, audit and execution | Interface compatibility |
| `cases.py` | Materialized example/source case dictionaries | Parameter provenance |
| `schedules.py` | Bed/step topology and connection semantics | Cycle physics |
| `core.py` | State layout, equilibrium, LDF, balances, valves and network RHS | Highest scientific risk |
| `numerics.py` | Jacobian and integration utilities | Stability and performance |
| `vectorized.py` | Vectorized physics paths | Scalar/vector agreement |
| `runner.py` | BDF solve, ledgers, CSS, checkpoint/restart and run artifacts | Reproducibility |
| `qualification.py` | Input and trajectory domain screens | Evidence boundary |
| `release_filter.py` | Purpose-specific default-deny result selection | Downstream safety |
| `mean_stream_bridge.py` | Cycle-average FEED/PRODUCT/TAIL contract | Process coupling |
| `phase_aligned_state.py` | Phase-consistent stream state construction | Thermodynamic semantics |
| `dwsim_payload_preflight.py` | Reject incomplete or scientifically closed payloads | DWSIM boundary |
| `peng_robinson.py` | Optional real-gas closure | Property-method validity |
| `guderian_2021_n2.py` | Source-isolated N2 PSA adapter | Source-specific only |
| `dantas_2011_co2_n2_adapter.py` | Source-isolated CO2/N2 adapter | Source-specific only |
| `tamof1_2025_adapter.py` | Source-isolated TAMOF-1 adapter | Source-specific only |
| `streb_2021_*.py` | Source-specific pressure, equilibrium and timeline helpers | Source-specific only |
| `batch.py` | Immutable independent-case subprocess dispatch | Campaign isolation |

## State and result identity

`runner.py` defines the engine files that participate in the restart identity.
Changing a state-transition, schedule or qualification module changes the
engine hash and makes an old checkpoint ineligible for silent continuation.
Audit/report-only changes should not invalidate a long physical trajectory,
but the full source snapshot is still archived with the run.

Important artifacts:

| Artifact | Meaning |
|---|---|
| `input.json` | Frozen materialized case |
| `schedule.json` | Executed topology and timing |
| `engineering_input_screen.json` | Missing source/design evidence before solving |
| `source_code.zip` | Source snapshot for provenance |
| `checkpoint.json` | Committed resume metadata and state hash |
| `checkpoint.state.0/1.npz` | Alternating process-state generations |
| `active_step.json` | Non-terminal progress heartbeat |
| cycle/step ledgers | Integrated boundary and balance evidence |
| `result.json` | Terminal or deliberately paused run state |
| `failure.json` | Structured caught failure |
| `failure_history/` | Preserved earlier failures before a safe retry |

## Test map

| Test file | Main contract |
|---|---|
| `test_peng_robinson.py` | Real-gas calculations and derivatives |
| `test_vectorized.py` | Scalar/vectorized physics agreement |
| `test_schedule_topology_contract.py` | Bed/step and connection semantics |
| `test_qualification.py` | Fail-closed engineering/domain screens |
| `test_batch.py` | Independent-case dispatch |
| `test_batch_hardening.py` | Path, hash and failure isolation |
| `test_release_filter.py` | Explicit boolean release policy |
| `test_mean_stream_bridge.py` | Cycle-average stream accounting |
| `test_phase_aligned_state.py` | Phase-aligned state construction |
| `test_dwsim_payload_preflight.py` | Rejection of incomplete/unreleased transfers |
| `test_guderian_2021_*.py` | Source-specific N2 adapter and schedule |

When changing a module, run the whole suite. Scientific interactions make a
narrow unit-test pass insufficient for release.

## Public evidence files

- `validation/validation_summary.json` is the compact machine-readable metric
  record.
- `docs/research-state.json` is the project-level gate and next-action record.
- `docs/validation.md` explains the evidence in prose.
- `docs/known-issues.md` preserves failed and blocked approaches.
- `docs/references.md` maps source families to their permitted use.

The original large calculation directories are not public. A public claim must
therefore be supportable by the committed derived metrics, tests and source
locators. If a stronger claim needs a non-public artifact, state that
dependency and keep the release gate closed.

## Safe extension pattern

For a new adsorbent or paper:

1. add a dedicated source adapter rather than editing an unrelated one;
2. add the source contract and unit conversions;
3. add a small parse/topology test;
4. create calibration and holdout cases before scoring;
5. run numerical qualification;
6. produce a machine-readable decision with all downstream flags false by
   default;
7. update the validation summary and research state only after review.

