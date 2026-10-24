# Input and output reference

## Case schema

The top-level input schema is `dynamic-psa-case/v1`.

```json
{
  "schema": "dynamic-psa-case/v1",
  "case_id": "example",
  "species": [],
  "bed": {},
  "feed": {},
  "high_pressure_bar": 15.65,
  "low_pressure_bar": 1.1,
  "product_pressure_bar": 15.0,
  "initial_pressure_bar": 1.1,
  "purge_to_feed_ratio": 0.05,
  "boundary_mode": "compressible_valve",
  "valves": {},
  "schedule": {},
  "solver": {},
  "provenance": {},
  "claim": ""
}
```

### Species

Each species supplies `name`, `mw_kg_mol` and `cp_j_mol_k`. Peng-Robinson
cases additionally require critical temperature, critical pressure and
acentric factor for every component.

### Bed

Important fields include:

- `n_cells`, `length_m`, `diameter_m`;
- interparticle and optional pellet void fraction;
- adsorbent density and particle diameter;
- solid, bed and wall thermal properties;
- bed-wall and ambient heat-transfer coefficients;
- `kinetics`: `loading_ldf` or `partial_pressure_ldf`;
- `energy_mode`;
- one or more `materials`, ordered by increasing `end_fraction`.

### Feed

`mole_fractions`, `temperature_k`, `pressure_bar` and `molar_flow_mol_s` define
the feed boundary. Mole fractions must match the species order.

### Schedule

The schedule name selects a topology plugin. Timing fields are specific to the
plugin. For example, `skarstrom_4step` uses pressurization and adsorption
times, whereas `equalized_8step` uses a slot time.

### Solver

The solver section controls relative tolerance, maximum step, maximum and
minimum cycles, CSS tolerance, consecutive CSS cycles, balance tolerances and
optional trace recording.

## Output directory

| File | Meaning |
|---|---|
| `input.json` | Frozen case actually consumed |
| `schedule.json` | Materialized bed/slot schedule |
| `engineering_input_screen.json` | Pre-run applicability and evidence screen |
| `checkpoint.json` | Restart metadata |
| `checkpoint.state.*.npz` | Alternating numerical state checkpoints |
| `active_step.json` | Live physical-slot progress |
| `cycles.csv` | One row per completed cycle |
| `steps.csv` | One row per completed physical slot |
| `last_cycle_trace.csv` | Retained samples from the latest cycle |
| `trajectory_domain_screen.json` | Post-run material-domain screen |
| `result.json` | Terminal result and release flags |
| `failure.json` | Structured failure; never a valid result row |
| `source_code.zip` | Deterministic source snapshot |

## Important result fields

- `status`;
- `state_css_converged`;
- `number_of_beds` and `number_of_cycle_steps`;
- `final_cycle.reporting_purity` and `reporting_recovery`;
- `final_cycle.mass_balance_relative`;
- `final_cycle.energy_balance_relative`;
- component feed/product/tail ledgers;
- enthalpy, ambient-heat and vacuum-work ledgers when available;
- release and training eligibility flags.

Do not use a wildcard to collect every `result.json` as training data. Only an
explicit allow-list produced by the release audit is admissible.
