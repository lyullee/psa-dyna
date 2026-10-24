# Physics and numerical formulation

## 1. Modeling basis

PSA-Dyna represents each packed bed as a one-dimensional axial domain. The
domain is divided into finite volumes. The solver evolves conservative cell
inventories rather than mole fractions. Composition is reconstructed from
non-negative component inventories.

This choice makes external feed, product, tail gas and inter-bed pressure
equalization use the same component balance. Internal equalization is never
counted as fresh feed.

## 2. Component conservation

For component `i` in cell `j`, the semi-discrete balance has the form

```text
d n(i,j)/dt = F(i,j-1/2) - F(i,j+1/2) - m_s(j) d q(i,j)/dt
```

where `n` is gas inventory, `F` is component molar flux, `m_s` is adsorbent
mass and `q` is adsorbed loading. The face flux contains an upwind advective
term and, when enabled, an axial dispersion term.

## 3. Gas pressure and equation of state

The default gas closure is ideal gas:

```text
P V = n R T
```

An optional Peng-Robinson path calculates density from pressure, temperature
and composition. When explicitly enabled, departure internal energy and stream
enthalpy are included consistently. Critical properties and binary interaction
parameters must be supplied by the case; the engine does not invent them.

## 4. Packed-bed momentum closure

Internal face velocity is obtained from the Ergun relation. The code solves the
quadratic pressure-gradient form without subtractive cancellation in the Darcy
limit. Local voidage, particle diameter, viscosity and gas density determine
the linear and inertial resistance.

The face molar flux uses the donor-cell composition. Reverse internal flow is
therefore handled without assigning a negative composition.

## 5. Adsorption equilibrium

The standard multicomponent equilibrium option is competitive dual-site
Langmuir. For site `s`:

```text
q*(i,s) = qsat(i,s,T) b(i,s,T) p_i /
          [1 + sum_k b(k,s,T) p_k]
```

with temperature-dependent affinity and, where specified, temperature-
dependent capacity. Analytical pressure and temperature derivatives are used
in the coupled pellet/energy calculation.

Source-local adapters implement additional published forms, including the
Guderian Sips/IAST formulation. These parameter sets remain isolated; values
from different adsorbents or gas systems must not be mixed.

## 6. Adsorption kinetics

Two distinct kinetic state formulations are supported.

### Loading LDF

```text
dq_i/dt = k_i [q*_i(P,T,y) - q_i]
```

This form stores loading directly and is used by the layered AC/5A process
screen and source-specific loading models.

### Pellet partial-pressure LDF

The Wang/RB3 apparatus model stores a pellet pressure-like state. Gas-pellet
exchange is driven by the difference between gas partial pressure and the
pellet state. Equilibrium derivatives couple pellet accumulation and
temperature. The two kinetic forms are not aliases and are not interchanged.

## 7. Energy balances

The bed energy ledger includes, as selected by the case:

- gas sensible internal energy;
- optional Peng-Robinson departure energy;
- solid thermal capacity;
- adsorption energy;
- axial bed conduction;
- bed-to-wall heat transfer;
- axial wall conduction;
- wall-to-ambient heat transfer;
- enthalpy carried by every boundary and internal face flow.

For a dynamic wall, the wall temperature is an independent state. A
fixed-ambient wall option is also available for source-specific experiments.

The cell energy equation is solved together with adsorption states. Adsorption
heat is therefore coupled to equilibrium, kinetics and gas inventory rather
than applied as an external temperature correction.

## 8. Valves and external boundaries

The general valve model is a one-way compressible orifice law. It calculates
subcritical or choked flow from upstream pressure, downstream pressure,
temperature, mixture molecular weight, heat-capacity ratio and `CdA`.

Near zero pressure difference, a small regularization removes the infinite
slope of the inviscid square-root relation. Reverse pressure difference closes
the one-way valve exactly. Source-specific controllers are kept as separately
named boundary modes.

## 9. Multi-bed cycle scheduling

The schedule defines a matrix of bed actions and physical time slots. Bed
count and step count are independent metadata and are verified against the
matrix before integration.

Included engineering schedules are:

- two-bed, four-step Skarstrom cycle;
- four-bed, eight-step cycle with two pressure-equalization levels.

Source-specific schedules are labeled separately when published valve timing
is incomplete or inferred.

## 10. Numerical integration and CSS

The finite-volume equations form a stiff ODE system integrated with SciPy.
The runner records a deterministic source snapshot, frozen input hash,
alternating checkpoints and cycle ledgers.

Ordinary cyclic steady state requires three consecutive cycles satisfying all
of the following default limits:

- full-state scaled change <= `1e-4`;
- pressure and temperature scaled change <= `1e-4`;
- purity and recovery change <= `1e-4`;
- component and energy residuals <= `1e-5`.

Stable purity alone is not CSS.

## 11. Physical safeguards

The solver and post-processors check:

- finite, non-negative inventories and loadings;
- positive temperature and pressure;
- material fit-domain membership;
- component and energy conservation;
- schedule topology consistency;
- input and source-code hashes;
- restart identity;
- downstream release eligibility.

Failure of a safeguard creates an explicit failure or closed-release record;
it is not silently converted into a valid data row.
