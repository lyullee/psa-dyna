# Known limitations

1. The current public experimental pass is a single-bed, figure-level
   breakthrough comparison. It is not a cyclic industrial PSA validation.
2. Public publisher figures do not provide the same information as raw,
   synchronized pressure, temperature, composition and valve histories.
3. Several cycle timings and valve areas in engineering schedules are design
   assumptions rather than equipment-vendor data.
4. The 1D model cannot resolve radial temperature, channeling, maldistribution,
   attrition, pellet cracking or detailed distributor hydrodynamics.
5. Water, sulfur and other trace contaminants require qualified adsorption and
   regeneration models. A dry-gas result must not be treated as a wet-biogas
   guard-bed design.
6. The default gas model is ideal. Peng-Robinson is optional but still requires
   validated component and interaction data.
7. Heat-transfer coefficients are case inputs. Scaling a laboratory `k/r`
   estimate to a large vessel does not establish plant-scale thermal behavior.
8. Numerical CSS does not establish material applicability or mechanical
   equipment feasibility.
9. The DWSIM bridge is a guarded JSON contract, not a released automatic
   flowsheet writer.
10. The software is a research prototype and is not a pressure-vessel,
    functional-safety or operating-control design tool.
