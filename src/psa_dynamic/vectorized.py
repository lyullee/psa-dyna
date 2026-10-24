"""Batched finite-volume algebra, with the scalar core retained as a reference.

Vectorization changes neither the constitutive laws nor the upwind stencil.
Arrays retain complex perturbations for the active-branch Jacobian. Materials
are batched separately, so a layer may have a different number of sites.
"""
import numpy as np
from .qualification import check_capacity

R = 8.314462618
BAR = 1.e5


def equilibrium_batch(pressure, temperature, material):
    """Return q[cell,species,site], dq/dp[cell,site,i,j] and dq/dT.

    This is the same competitive Langmuir algebra as core.equilibrium, not a
    surrogate. The real branch controls positivity; imaginary parts survive.
    """
    p = np.where(np.real(pressure) >= 0., pressure, 0.)
    t = np.asarray(temperature)[:, None, None]
    if material.get("model") in {
        "streb_sips_wilson_rast_candidate",
        "guderian_2021_sips_iast",
    }:
        # RAST and IAST both require a nonlinear mixture root per cell.  Reuse
        # the scalar source-local mathematics exactly; substituting a
        # competitive Langmuir vectorization would change the physical model.
        # These loading-state branches do not consume equilibrium derivatives.
        from .core import equilibrium
        q = np.stack([equilibrium(p[i], temperature[i], material)[0]
                      for i in range(len(p))], axis=0)
        return q, None, None
    if material.get("model") == "extended_langmuir_3":
        ip1 = np.asarray(material["ip1_kmol_kg"], float)[None, :, :] * 1000.0
        ip2 = np.asarray(material["ip2_kmol_kg_k"], float)[None, :, :] * 1000.0
        ip3 = np.asarray(material["ip3_per_bar"], float)[None, :, :]
        ip4 = np.asarray(material["ip4_k"], float)[None, :, :]
        if any(arr.shape != ip1.shape for arr in (ip2, ip3, ip4)):
            raise ValueError("Inconsistent Extended Langmuir-3 parameter dimensions")
        qmax = ip1 - ip2 * t
        if np.any(np.real(qmax) <= 0.):
            raise ValueError("Extended Langmuir-3 qmax is non-positive")
        p_safe = np.where(np.real(p) > 0., p, 1e-30)
        u = ip3 * np.exp(ip4 / t) * p_safe[:, :, None]
        den = 1. + np.sum(u, axis=1, keepdims=True)
        q = qmax * u / den
        b = ip3 * np.exp(ip4 / t)
        du_dp = b
        du_dt = u * (-ip4 / t**2)
        dden_dt = np.sum(du_dt, axis=1, keepdims=True)
        qt = -ip2 * u / den + qmax * (du_dt / den - u * dden_dt / den**2)
        diagonal = (qmax * du_dp / den).transpose(0, 2, 1)
        left = (qmax * u / den**2).transpose(0, 2, 1)
        right = du_dp.transpose(0, 2, 1)
        jac = diagonal[:, :, :, None] * np.eye(p.shape[1]) - left[:, :, :, None] * right[:, :, None, :]
        return q, jac, qt
    if material.get("model") == "extended_langmuir_freundlich":
        # Same source-labelled extended Langmuir--Freundlich law as the
        # scalar implementation.  The explicit atm/bar scale is part of the
        # material record because the published Jee table uses atm.
        k1 = np.asarray(material["lf_k1_mol_kg"], float)[None, :, :]
        k2 = np.asarray(material["lf_k2_mol_kg_k"], float)[None, :, :]
        k3 = np.asarray(material["lf_k3_per_atm_power"], float)[None, :, :]
        k4 = np.asarray(material["lf_k4_k"], float)[None, :, :]
        k5 = np.asarray(material["lf_k5"], float)[None, :, :]
        k6 = np.asarray(material["lf_k6_k"], float)[None, :, :]
        if any(arr.shape != k1.shape for arr in (k2, k3, k4, k5, k6)):
            raise ValueError("Inconsistent extended Langmuir-Freundlich parameter dimensions")
        scale = float(material.get("lf_pressure_scale_bar", 1.01325))
        qmax = k1 + k2 * t
        if np.any(np.real(qmax) <= 0.):
            raise ValueError("Extended Langmuir-Freundlich qmax is non-positive")
        p_atm = p[:, :, None] / scale
        p_safe = np.where(np.real(p_atm) > 0., p_atm, 1e-30)
        b = k3 * np.exp(k4 / t)
        exponent = k5 + k6 / t
        u = b * p_safe ** exponent
        den = 1. + np.sum(u, axis=1, keepdims=True)
        q = qmax * u / den
        du_dp = exponent * b * p_safe ** (exponent - 1.) / scale
        du_dt = u * (-k4 / t**2 - k6 * np.log(p_safe) / t**2)
        dden_dt = np.sum(du_dt, axis=1, keepdims=True)
        qt = k2 * u / den + qmax * (du_dt / den - u * dden_dt / den**2)
        diagonal = (qmax * du_dp / den).transpose(0, 2, 1)
        left = (qmax * u / den**2).transpose(0, 2, 1)
        right = du_dp.transpose(0, 2, 1)
        jac = diagonal[:, :, :, None] * np.eye(p.shape[1]) - left[:, :, :, None] * right[:, :, None, :]
        return q, jac, qt
    if material.get("model") == "dantas_toth_source":
        # Source-local independent pure-component Toth branch. Dantas' binary
        # runs use separate CO2 and N2 fits plus an explicit inert He carrier;
        # do not introduce a competitive Langmuir denominator here.
        qm = np.asarray(material["dantas_qm_mol_kg"], float)[None, :, :]
        exponent = np.asarray(material["dantas_n"], float)[None, :, :]
        k0 = np.asarray(material["dantas_k0_per_bar"], float)[None, :, :]
        heat = np.asarray(material["dantas_heat_j_mol"], float)[None, :, :]
        if qm.shape[1:] != (p.shape[1], 1) or any(value.shape != qm.shape for value in (exponent, k0, heat)):
            raise ValueError("Inconsistent Dantas Toth parameter dimensions")
        if np.any(np.real(qm) < 0.0) or np.any(np.real(exponent) <= 0.0):
            raise ValueError("Dantas Toth parameters require non-negative capacity and positive exponent")
        p_safe = np.where(np.real(p) > 0.0, p, 1e-30)
        t_cells = t[:, 0, 0][:, None, None]
        # Dantas Table 6 reports positive -DeltaH; Eq. (18) therefore gives
        # K=K0*exp(+|DeltaH|/(R*T)).
        k = k0 * np.exp(heat / (R * t_cells))
        u = k * p_safe[:, :, None]
        un = u ** exponent
        base = 1.0 + un
        q = qm * u / base ** (1.0 / exponent)
        dq_dp = qm * k * base ** (-1.0 / exponent - 1.0)
        du_dt = -u * heat / (R * t_cells**2)
        qt = qm * base ** (-1.0 / exponent - 1.0) * du_dt
        jac = np.zeros((p.shape[0], p.shape[1], p.shape[1]), dtype=np.result_type(p, temperature, float))
        diagonal = dq_dp[:, :, 0]
        indices = np.arange(p.shape[1])
        jac[:, indices, indices] = diagonal
        return q, jac, qt
    if material.get("model") == "tian_toth_source":
        # Tian 2018 pure-component Toth namespace.  No competitive mixture
        # denominator is introduced because the source does not publish one.
        # n(T)=t0+c/T is retained so the thermal derivative is source-faithful.
        qm = np.asarray(material["tian_qm_mol_kg"], float)[None, :, :]
        t0 = np.asarray(material["tian_t0"], float)[None, :, :]
        c_k = np.asarray(material["tian_c_k"], float)[None, :, :]
        k0 = np.asarray(material["tian_k0_per_bar"], float)[None, :, :]
        heat = np.asarray(material["tian_heat_j_mol"], float)[None, :, :]
        if qm.shape[1:] != (p.shape[1], 1) or any(value.shape != qm.shape for value in (t0, c_k, k0, heat)):
            raise ValueError("Inconsistent Tian Toth parameter dimensions")
        if np.any(np.real(qm) < 0.0):
            raise ValueError("Tian Toth parameters require non-negative capacity")
        p_safe = np.where(np.real(p) > 0.0, p, 1e-30)
        t_cells = t[:, 0, 0][:, None, None]
        exponent = t0 + c_k / t_cells
        if np.any(np.real(exponent) <= 0.0):
            raise ValueError("Tian Toth temperature-dependent exponent must be positive")
        k = k0 * np.exp(heat / (R * t_cells))
        u = k * p_safe[:, :, None]
        u_safe = np.where(np.real(u) > 0.0, u, 1e-300)
        un = u_safe ** exponent
        base = 1.0 + un
        q = qm * u / base ** (1.0 / exponent)
        dq_dp = qm * k * base ** (-1.0 / exponent - 1.0)
        dlnu_dt = -heat / (R * t_cells**2)
        dn_dt = -c_k / t_cells**2
        dlnbase_dt = (un / base) * (dn_dt * np.log(u_safe) + exponent * dlnu_dt)
        dlnq_dt = dlnu_dt + (dn_dt / exponent**2) * np.log(base) - dlnbase_dt / exponent
        qt = q * dlnq_dt
        jac = np.zeros((p.shape[0], p.shape[1], p.shape[1]), dtype=np.result_type(p, temperature, float))
        diagonal = dq_dp[:, :, 0]
        indices = np.arange(p.shape[1])
        jac[:, indices, indices] = diagonal
        return q, jac, qt
    if material.get("model") == "tamof1_dual_site_langmuir_freundlich_303k":
        # TAMOF-1 source-local dual-site Langmuir--Freundlich equilibrium.
        # Supplementary Table 16 is a 303 K equilibrium table, not the
        # legacy a1/a2/b0 correlation.  Keep this branch before the legacy
        # fallback so the source parameters can never be silently mixed.
        qmax = np.asarray(material["tamof_qm_mol_kg"], dtype=float)
        affinity = np.asarray(material["tamof_k_per_bar_power"], dtype=float)
        exponent = np.asarray(material["tamof_n"], dtype=float)
        if qmax.shape != (p.shape[1], 2) or affinity.shape != qmax.shape or exponent.shape != qmax.shape:
            raise ValueError("TAMOF-1 dual-site LF requires species x two-site arrays")
        reference = float(material.get("tamof_reference_temperature_k", 303.0))
        tolerance = float(material.get("tamof_temperature_tolerance_k", 1.0e-9))
        if np.any(np.abs(np.real(t[:, 0, 0]) - reference) > tolerance):
            raise ValueError("TAMOF-1 303 K equilibrium is outside its published temperature domain")
        if np.any(qmax < 0.0) or np.any(affinity <= 0.0) or np.any(exponent <= 0.0):
            raise ValueError("TAMOF-1 dual-site LF parameters must be non-negative/positive")
        p_safe = np.where(np.real(p) > 0.0, p, 1.0e-30)
        u = (affinity[None, :, :] * p_safe[:, :, None]) ** exponent[None, :, :]
        den = 1.0 + np.sum(u, axis=1, keepdims=True)
        q = qmax[None, :, :] * u / den
        jac = np.empty((p.shape[0], 2, p.shape[1], p.shape[1]), dtype=np.result_type(p, temperature, float))
        qt = np.zeros_like(q)
        indices = np.arange(p.shape[1])
        for site in range(2):
            du_dp = exponent[None, :, site] * u[:, :, site] / p_safe
            denominator = den[:, 0, site][:, None]
            diagonal = qmax[None, :, site] * du_dp / denominator
            left = qmax[None, :, site] * u[:, :, site] / denominator ** 2
            jac[:, site] = -left[:, :, None] * du_dp[:, None, :]
            jac[:, site, indices, indices] += diagonal
        return q, jac, qt
    a1 = np.asarray(material["a1"], float)
    a2 = np.asarray(material.get("a2", np.zeros_like(a1)), float)
    energy = np.asarray(material.get("affinity_energy_j_mol", material["binding_energy_j_mol"]), float)
    b = np.asarray(material["b0_per_bar"], float) * np.exp(energy / (R*t))
    cap = a1 - a2*t
    check_capacity(material,temperature,cap)
    den = 1. + np.sum(b*p[:, :, None], axis=1, keepdims=True)
    q = cap*b*p[:, :, None]/den
    db = -b*energy/(R*t*t)
    dden = np.sum(db*p[:, :, None], axis=1, keepdims=True)
    qt = -a2*b*p[:, :, None]/den + cap*p[:, :, None]*(db/den-b*dden/den**2)
    # Transpose the site dimension before forming the species-by-species
    # outer product. Each derivative includes competition at the same site.
    diagonal = (cap*b/den).transpose(0, 2, 1)
    left = (cap*b*p[:, :, None]/den**2).transpose(0, 2, 1)
    right = b.transpose(0, 2, 1)
    jac = diagonal[:, :, :, None]*np.eye(p.shape[1]) - left[:, :, :, None]*right[:, :, None, :]
    return q, jac, qt


def internal_flux(bed, state):
    """Compute all interior Ergun/upwind faces in one NumPy batch."""
    gas, _, t, _, p, y = bed.fields(state)
    dtype = np.result_type(state, float)
    flux = np.zeros((bed.n+1, bed.ns), dtype=dtype)
    heat = np.zeros(bed.n+1, dtype=dtype)
    eps = 0.5 * (bed.void_fraction[:-1] + bed.void_fraction[1:])
    particle = 0.5 * (bed.particle_diameter[:-1] + bed.particle_diameter[1:])
    linear = 150*bed.spec.get("viscosity_pa_s", 1.2e-5)*(1-eps)**2/(eps**3*particle**2)
    delta = (p[:-1]-p[1:])*BAR
    pf = (p[:-1]+p[1:])*BAR/2
    tf = (t[:-1]+t[1:])/2
    yf = (y[:-1]+y[1:])/2
    if bed.equation_of_state == "peng_robinson":
        concentration = 0.5 * (gas[:-1].sum(axis=1) / bed.vgas[:-1]
                               + gas[1:].sum(axis=1) / bed.vgas[1:])
        rho = concentration * (yf @ bed.mw)
    else:
        rho = pf*(yf@bed.mw)/(R*tf)
    quadratic = 1.75*rho*(1-eps)/(eps**3*particle)
    sign = np.where(np.real(delta) >= 0., 1., -1.)
    gradient = sign*delta/bed.face_spacing
    velocity = 2*gradient/(linear+np.sqrt(linear**2+4*quadratic*gradient))
    if bed.equation_of_state == "peng_robinson":
        flow = sign * velocity * bed.area * concentration
    else:
        flow = sign*velocity*bed.area*pf/(R*tf)
    donor_y = np.where((sign > 0.)[:, None], y[:-1], y[1:])
    donor_t = np.where(sign > 0., t[:-1], t[1:])
    reference = 0. if bed.energy_mode == "internal_energy" else bed.spec["ambient_k"]
    flux[1:-1] = flow[:, None]*donor_y
    if bed.pr_departure_energy:
        donor_p = np.where(sign > 0.0, p[:-1], p[1:])
        donor_h = np.asarray([
            bed.enthalpy(donor_y[index], donor_t[index], donor_p[index])
            for index in range(len(flow))
        ])
        heat[1:-1] = flow * donor_h
    else:
        heat[1:-1] = flow*(donor_y@bed.cp)*(donor_t-reference)
    if bed.axial_dispersion_model == "guderian_2021_eq3":
        from .guderian_2021_n2 import axial_dispersion_o2_n2_m2_s
        dispersion = np.stack([
            axial_dispersion_o2_n2_m2_s(
                velocity[index], pf[index]/BAR, tf[index], eps[index],
                bed.dispersion_particle_radius_m,
            )
            for index in range(len(velocity))
        ])
        concentration = 0.5 * (gas[:-1].sum(axis=1) / bed.vgas[:-1]
                               + gas[1:].sum(axis=1) / bed.vgas[1:])
        diffusive = (-dispersion * bed.area * concentration[:, None]
                     * (y[1:] - y[:-1]) / bed.face_spacing[:, None])
        flux[1:-1] += diffusive
        yf = 0.5 * (y[:-1] + y[1:])
        tf = 0.5 * (t[:-1] + t[1:])
        if bed.pr_departure_energy:
            diffusive_h = np.asarray([
                bed.enthalpy(yf[index], tf[index], pf[index] / BAR)
                for index in range(len(flow))
            ])
            heat[1:-1] += diffusive.sum(axis=1) * diffusive_h
        else:
            heat[1:-1] += diffusive.sum(axis=1) * np.einsum("ij,j->i", yf, bed.cp) * (tf-reference)
    elif bed.axial_dispersion_m2_s > 0:
        concentration = 0.5 * (gas[:-1].sum(axis=1) / bed.vgas[:-1]
                               + gas[1:].sum(axis=1) / bed.vgas[1:])
        diffusive = (-bed.axial_dispersion_m2_s * bed.area * concentration[:, None]
                     * (y[1:] - y[:-1]) / bed.face_spacing[:, None])
        flux[1:-1] += diffusive
        yf = 0.5 * (y[:-1] + y[1:])
        tf = 0.5 * (t[:-1] + t[1:])
        if bed.pr_departure_energy:
            diffusive_h = np.asarray([
                bed.enthalpy(yf[index], tf[index], pf[index] / BAR)
                for index in range(len(flow))
            ])
            heat[1:-1] += diffusive.sum(axis=1) * diffusive_h
        else:
            heat[1:-1] += diffusive.sum(axis=1) * np.einsum("ij,j->i", yf, bed.cp) * (tf-reference)
    return flux, heat


def rates(bed, state, flux, enthalpy_flux, h_bed_wall_override=None):
    """Batch cell balances and small pellet/temperature linear systems."""
    gas, ads, t, tw, p, y = bed.fields(state)
    flow_rate = flux[:-1]-flux[1:]
    heat = enthalpy_flux[:-1]-enthalpy_flux[1:]
    dtype = np.result_type(state, float)
    qbed = np.zeros(bed.n, dtype=dtype)
    qwall = np.zeros(bed.n, dtype=dtype)
    bed_k_face = 0.5 * (bed.bed_k[:-1] + bed.bed_k[1:])
    conduction_factor_face = 0.5 * (bed.conduction_factor[:-1] + bed.conduction_factor[1:])
    qb = bed_k_face*conduction_factor_face*bed.area/bed.face_spacing*np.diff(t)
    qw = bed.spec.get("wall_k_w_m_k", 10.9)*bed.wall_cross_area/bed.face_spacing*np.diff(tw)
    qbed[:-1] += qb
    qbed[1:] -= qb
    qwall[:-1] += qw
    qwall[1:] -= qw
    # A source-labelled cycle may specify a different internal gas-wall
    # coefficient for each operation.  Passing it per RHS evaluation keeps
    # the shared Bed immutable across beds and solver/Jacobian calls.
    h_bed_wall = (bed.spec.get("h_bed_wall_w_m2_k", 0.) if h_bed_wall_override is None
                  else h_bed_wall_override)
    if bed.wall_model == "fixed_ambient":
        wall_exchange = h_bed_wall*bed.wall_area*(bed.spec["ambient_k"]-t)
        ambient = wall_exchange.copy()
    else:
        wall_exchange = h_bed_wall*bed.wall_area*(tw-t)
        ambient = bed.spec.get("h_ambient_w_m2_k", 0.)*bed.outer_area*(bed.spec["ambient_k"]-tw)
    dx = np.zeros((bed.n, bed.nv), dtype=dtype)
    if bed.wall_model == "fixed_ambient":
        dx[:, -1] = 0.0
    else:
        dx[:, -1] = (-wall_exchange+qwall+ambient)/bed.cwall
    # Mesh construction makes each material contiguous. No run-time search
    # or material interpolation is introduced at a layer boundary.
    for cells, material in bed.material_batches:
        model = material.get("model", "dual_site_langmuir")
        tj, aj, gj = t[cells], ads[cells], gas[cells]
        mass = bed.mass[cells]
        if bed.mode == "loading_ldf":
            if model == "guderian_2021_sips_iast":
                q, _, _ = equilibrium_batch(p[cells, None]*y[cells], tj, material)
                from .guderian_2021_n2 import (
                    SPECIES_ORDER,
                    adsorption_energy_temperature_derivative_j_kg_k,
                    heat_of_adsorption_magnitude_j_mol,
                    qdf_rate_mol_kg_s,
                )

                target = q.sum(axis=2)
                da = np.stack([
                    qdf_rate_mol_kg_s(aj[index], target[index], tj[index])
                    for index in range(len(tj))
                ])
                heat_weights = np.asarray([
                    [heat_of_adsorption_magnitude_j_mol(
                        name, aj[cell_index, species_index], tj[cell_index]
                    ) for species_index, name in enumerate(SPECIES_ORDER)]
                    for cell_index in range(len(tj))
                ])
                energy_temperature_derivative = np.asarray([
                    sum(adsorption_energy_temperature_derivative_j_kg_k(
                        name, aj[cell_index, species_index], tj[cell_index]
                    ) for species_index, name in enumerate(SPECIES_ORDER))
                    for cell_index in range(len(tj))
                ])
            else:
                weights = np.asarray(material["binding_energy_j_mol"])
                ldf = np.asarray(material["ldf_per_s"])
                heat_weights = np.broadcast_to(weights[:, 0], aj.shape)
                energy_temperature_derivative = np.zeros(len(tj), dtype=dtype)
                # A zero-LDF material is an explicit non-adsorbing/inert
                # transport reference.  Do not evaluate an equilibrium model
                # that cannot influence the derivative: during stiff solver
                # trial states ``0 * inf`` otherwise creates NaNs even though
                # the mathematically exact loading rate is zero.
                if np.all(ldf == 0.0):
                    da = np.zeros_like(aj)
                else:
                    q, _, _ = equilibrium_batch(p[cells, None]*y[cells], tj, material)
                    da = ldf*(q.sum(axis=2)-aj)
            dn = flow_rate[cells]-mass[:, None]*da
            thermal_mass = bed.thermal_mass[cells]
            rhs = (heat[cells]+qbed[cells]+wall_exchange[cells]
                   + thermal_mass*np.sum(da*heat_weights, axis=1))
            capacity = (bed.thermal_cs[cells]
                        - thermal_mass*energy_temperature_derivative)
            if bed.energy_mode == "internal_energy":
                if bed.pr_departure_energy:
                    component_energy = np.empty_like(gj)
                    gas_capacity = np.empty(len(tj), dtype=dtype)
                    cell_indices = np.arange(bed.n)[cells]
                    for local_index, cell_index in enumerate(cell_indices):
                        component_energy[local_index], gas_capacity[local_index] = (
                            bed.gas_internal_energy_derivatives(
                                gj[local_index], tj[local_index], int(cell_index)
                            )
                        )
                    rhs -= np.sum(dn * component_energy, axis=1)
                    capacity = capacity + gas_capacity
                else:
                    rhs -= (dn@bed.cv)*tj
                    capacity = capacity+gj@bed.cv
            if np.any(np.real(capacity) <= 0.0):
                raise FloatingPointError("Guderian/source thermal capacity became non-positive")
            dt = rhs/capacity
        else:
            weights = np.asarray(material["binding_energy_j_mol"])
            ldf = np.asarray(material["ldf_per_s"])
            q, jac_site, qt_site = equilibrium_batch(aj, tj, material)
            jac = jac_site.sum(axis=1)
            qt = qt_site.sum(axis=2)
            exchange = bed.vparticle[cells, None]*ldf*(p[cells, None]*y[cells]-aj)*BAR/(R*tj[:, None])
            dn = flow_rate[cells]-exchange
            alpha = bed.alpha[cells]
            matrix = np.zeros((len(tj), bed.ns+1, bed.ns+1), dtype=dtype)
            matrix[:, :bed.ns, :bed.ns] = (alpha/tj)[:, None, None]*np.eye(bed.ns)+mass[:, None, None]*jac
            matrix[:, :bed.ns, -1] = -alpha[:, None]*aj/tj[:, None]**2+mass[:, None]*qt
            binding_p = np.einsum("is,nsij->nj", weights, jac_site)
            binding_t = np.sum(weights*qt_site, axis=(1, 2))
            thermal_mass = bed.thermal_mass[cells]
            matrix[:, -1, :bed.ns] = -thermal_mass[:, None]*binding_p
            matrix[:, -1, -1] = bed.thermal_cs[cells]-thermal_mass*binding_t
            rhs_heat = heat[cells]+qbed[cells]+wall_exchange[cells]
            if bed.energy_mode == "internal_energy":
                matrix[:, -1, :bed.ns] += alpha[:, None]*bed.cv
                matrix[:, -1, -1] += gj@bed.cv
                rhs_heat -= (dn@bed.cv)*tj
            rhs = np.column_stack([exchange, rhs_heat])
            # Explicit trailing RHS axis works with NumPy 1.x and 2.x alike.
            solved = np.linalg.solve(matrix, rhs[:, :, None])[:, :, 0]
            da, dt = solved[:, :-1], solved[:, -1]
        dx[cells, :bed.ns] = dn
        dx[cells, bed.ns:2*bed.ns] = da
        dx[cells, -2] = dt
    return dx.ravel(), ambient.sum()
