# S-015 minimum applicability supplement

This supplement answers one narrow question: does the current S-015
adsorbent and thermal model remain inside experimentally supported conditions?
It does not change model parameters or relax any previously declared gate.

## Accepted operating envelope

The accepted-node envelope from the ordinary N512 branch is:

| Quantity | Minimum | Maximum |
|---|---:|---:|
| 5A-layer temperature | 293.15 K | 368.59 K |
| 5A-layer total pressure | 1.10 bar | 15.14 bar |

These are marginal extrema from 210,172 accepted cell-time samples. The
temperature and pressure extremes need not occur in the same cell at the same
time.

## Adsorbent evidence map

The closest direct four-component equilibrium source is Yavary et al.
(`10.1021/acs.jced.6b00312`). It covers H2/CO2/CH4/CO at 298.15, 308.15 and
318.15 K over 1--18 bar. This supports a **screening comparison** for the
lower-temperature part of the S-015 path. It does not establish equivalence to
the Ahn commercial 5A, and the S-015 feed differs from the reported nominal
mixture by a total-variation distance of 0.106648.

The 318.15--368.59 K portion remains a hot-mixture extrapolation. Chen et al.
(`10.1016/0009-2509(90)80179-I`) is a relevant hot/high-pressure
multicomponent 5A scope anchor, but the public source does not provide a
replayable exact four-component table or prove equivalence to the Ahn grade.
Pure-component sources can bracket individual temperature or partial-pressure
axes, but they do not validate competitive cyclic adsorption.

Consequently, no screened source jointly establishes:

- exact Ahn adsorbent grade and binder;
- the full S-015 H2/CO2/CH4/CO mixture;
- the complete 293.15--368.59 K trajectory;
- cyclic LDF kinetics; and
- the plant geometry.

## Thermal evidence map

The case uses `h_bed_wall = 0.32 W m-2 K-1`. That value is exactly
`k_eff / R = 0.24 / 0.75`; it is a geometry-screening heuristic, not an
independently measured bed-to-wall heat-transfer coefficient. The vessel has
`D/dp = 543.48` and uses a one-dimensional axial temperature model.

Packed-bed experiments and models treat effective radial conductivity and the
wall coefficient as separate, condition-dependent parameters. Chan and
Lawther measured an air/glass-sphere bed and correlated wall transfer against
Reynolds number; Freiwald and Paterson examined the reliability of separately
estimated packed-bed heat-transfer parameters; elevated-pressure work also
measured radial conductivity and wall transfer separately. These sources
support the modelling requirement, not transfer of their coefficients into
S-015.

The phase-aligned N256-to-N512 RMS differences are 1.3468 K for gas
temperature and 1.4829 K for wall temperature. Both exceed the frozen 1 K
limits. Plant-scale thermal qualification therefore remains false.

## Minimum conclusion

The literature-domain screen is complete for a software-model paper. S-015 may
be described as a numerically converged, literature-parameterized screening
case with an explicit applicability map. It must not be described as an
experimentally qualified industrial PSA design.

No additional long CSS run can remove the missing exact-material and
plant-thermal evidence. HPC should be reserved for a preregistered N1024 or
thermal-sensitivity study only if the claim is expanded to plant-design
qualification.

Machine-readable derived evidence is provided in
[`validation/s015_minimum_applicability.json`](../validation/s015_minimum_applicability.json).

