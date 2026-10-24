"""Finite-volume PSA physics in SI inventories, bar pressures and kelvin.

States are *inventories*, never mole fractions. Bed transport, inter-bed
transfers and every external stream share the same RHS and time integrator.
This avoids the historical error of calling an internal recycle a new feed.
The two kinetic laws are deliberately distinct: loading LDF for the Ahn
design material and pellet partial-pressure LDF for the Wang apparatus.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.sparse import lil_matrix
from . import vectorized
from .peng_robinson import (
    PengRobinsonComponents,
    molar_enthalpy_departure as pr_molar_enthalpy_departure,
    molar_density_from_pressure as pr_molar_density_from_pressure,
    pressure_from_molar_density as pr_pressure_from_molar_density,
    total_internal_energy_departure as pr_total_internal_energy_departure,
    total_internal_energy_departure_derivatives as pr_total_internal_energy_departure_derivatives,
)
from .qualification import check_capacity
from .schedules import validate_schedule_contract
from .streb_2021_pressure import exponential_pressure_bar

R = 8.314462618
BAR = 1.0e5


def _positive(value, floor=0.):
    """Select the active real branch, preserving complex-step derivatives."""
    return np.where(np.real(value)>=floor,value,floor)


def _minimum(*values):
    return values[int(np.argmin([np.real(value) for value in values]))]


def equilibrium(p: np.ndarray, temperature: float, material: dict):
    """Competitive (dual-site) Langmuir and analytic p/T derivatives.

    All arrays have species/site axes. qmax=a1-a2*T [mol/kg];
    b=b0*exp(E/(R*T)) [1/bar]. The empirical affinity temperature
    coefficient need not equal an independently tabulated heat of adsorption.
    """
    p = _positive(p)
    if material.get("model") == "guderian_2021_sips_iast":
        # The source-local IAST root is real-valued and the Guderian adapter
        # uses loading QDF, so pressure/temperature Jacobians are not consumed.
        # Reject complex-step calls rather than silently dropping derivatives.
        if np.iscomplexobj(p) or np.iscomplexobj(temperature):
            raise TypeError("Guderian Sips/IAST requires a finite-difference Jacobian")
        from .guderian_2021_n2 import iast_loading_mol_kg
        q = iast_loading_mol_kg(
            np.asarray(p, dtype=float), float(temperature),
            tuple(material.get("species_order", ())),
        )
        return q[:, None], None, None
    if material.get("model") == "streb_sips_wilson_rast_candidate":
        # A nested nonlinear RAST root supplies loading only.  The present
        # adapter is restricted to loading LDF, which does not consume the
        # pressure/temperature equilibrium derivatives.  The standalone
        # solver is real-valued; complex-step Jacobians are rejected at the
        # Network level and here rather than silently losing imaginary parts.
        if np.iscomplexobj(p) or np.iscomplexobj(temperature):
            raise TypeError("Streb Wilson/RAST candidate requires a real-valued finite-difference Jacobian")
        from .streb_2021_wilson import solve_wilson_rast_mixture

        pressure = np.asarray(p, dtype=float)
        total = float(pressure.sum())
        if not math.isfinite(total) or total <= 0:
            raise ValueError("Streb Wilson/RAST candidate requires positive total pressure")
        # Cell inventories impose a tiny positive floor on absent species.
        # Preserve exact zero in the nonlinear active set below a numerical
        # 1e-14 mole-fraction threshold; this is far below product impurities.
        pressure = np.where(pressure / total < 1e-14, 0.0, pressure)
        active_total = float(pressure.sum())
        result = solve_wilson_rast_mixture(
            material["species_order"], pressure / active_total, active_total,
            float(temperature), material["wilson_lambda_matrix"],
            parameter_provenance=material["parameter_provenance"],
        )
        q = np.asarray(result["component_loadings_mol_kg"], dtype=float)[:, None]
        return q, None, None
    if material.get("model") == "extended_langmuir_3":
        # Aspen Adsorption Extended Langmuir-3 source form:
        # qmax_i=(IP1_i-IP2_i*T), b_i=IP3_i*exp(IP4_i/T),
        # q_i=qmax_i*b_i*p_i/(1+sum_j b_j*p_j).
        # IP1/IP2 are stored in kmol/kg in the source table and converted to
        # the engine's mol/kg inventory basis explicitly.
        ip1 = np.asarray(material["ip1_kmol_kg"], float) * 1000.0
        ip2 = np.asarray(material["ip2_kmol_kg_k"], float) * 1000.0
        ip3 = np.asarray(material["ip3_per_bar"], float)
        ip4 = np.asarray(material["ip4_k"], float)
        if any(arr.shape != ip1.shape for arr in (ip2, ip3, ip4)):
            raise ValueError("Inconsistent Extended Langmuir-3 parameter dimensions")
        t = temperature
        qmax = ip1 - ip2 * t
        if np.any(np.real(qmax) <= 0.):
            raise ValueError("Extended Langmuir-3 qmax is non-positive")
        p_safe = np.where(np.real(p) > 0., p, 1e-30)
        b = ip3 * np.exp(ip4 / t)
        u = b * p_safe[:, None]
        den = 1. + np.sum(u, axis=0)
        q = qmax * u / den
        du_dp = b
        du_dt = u * (-ip4 / t**2)
        dden_dt = np.sum(du_dt, axis=0)
        qt = -ip2 * u / den + qmax * (du_dt / den - u * dden_dt / den**2)
        jac = np.empty((ip1.shape[1], len(p), len(p)), dtype=np.result_type(p, temperature, float))
        for site in range(ip1.shape[1]):
            diagonal = qmax[:, site] * du_dp[:, site] / den[site]
            left = qmax[:, site] * u[:, site] / den[site]**2
            jac[site] = np.diag(diagonal) - left[:, None] * du_dp[:, site][None, :]
        return q, jac, qt
    if material.get("model") == "extended_langmuir_freundlich":
        # Source-labelled Jee/Cecchini cases use the extended
        # Langmuir--Freundlich form.  Pressures in the published table are
        # in atm, while the engine state is in bar; the explicit scale keeps
        # that unit conversion visible instead of silently changing B.
        k1 = np.asarray(material["lf_k1_mol_kg"], float)
        k2 = np.asarray(material["lf_k2_mol_kg_k"], float)
        k3 = np.asarray(material["lf_k3_per_atm_power"], float)
        k4 = np.asarray(material["lf_k4_k"], float)
        k5 = np.asarray(material["lf_k5"], float)
        k6 = np.asarray(material["lf_k6_k"], float)
        scale = float(material.get("lf_pressure_scale_bar", 1.01325))
        if any(arr.shape != k1.shape for arr in (k2, k3, k4, k5, k6)):
            raise ValueError("Inconsistent extended Langmuir-Freundlich parameter dimensions")
        t = temperature
        qmax = k1 + k2 * t
        if np.any(np.real(qmax) <= 0.):
            raise ValueError("Extended Langmuir-Freundlich qmax is non-positive")
        p_atm = p / scale
        p_safe = np.where(np.real(p_atm) > 0., p_atm, 1e-30)
        b = k3 * np.exp(k4 / t)
        exponent = k5 + k6 / t
        u = b * p_safe[:, None] ** exponent
        den = 1. + np.sum(u, axis=0)
        q = qmax * u / den
        du_dp = exponent * b * p_safe[:, None] ** (exponent - 1.) / scale
        du_dt = u * (-k4 / t**2 - k6 * np.log(p_safe[:, None]) / t**2)
        dden_dt = np.sum(du_dt, axis=0)
        qt = k2 * u / den + qmax * (du_dt / den - u * dden_dt / den**2)
        jac = np.empty((k1.shape[1], len(p), len(p)), dtype=np.result_type(p, temperature, float))
        for site in range(k1.shape[1]):
            diagonal = qmax[:, site] * du_dp[:, site] / den[site]
            left = qmax[:, site] * u[:, site] / den[site]**2
            jac[site] = np.diag(diagonal) - left[:, None] * du_dp[:, site][None, :]
        return q, jac, qt
    if material.get("model") == "competitive_sips_source":
        # Source-labelled competitive Sips form used by the Schell et al.
        # pre-combustion CO2/H2 benchmark.  The source reports A in 1/Pa,
        # while the engine stores partial pressure in bar; the conversion is
        # explicit so this parameter set cannot silently be mixed with the
        # engine's bar-based Langmuir tables.
        a = np.asarray(material["sips_a_mol_kg"], float)
        b = np.asarray(material["sips_b_j_mol"], float)
        A = np.asarray(material["sips_A_per_pa"], float)
        B = np.asarray(material["sips_B_j_mol"], float)
        alpha = np.asarray(material["sips_alpha"], float)
        beta = np.asarray(material["sips_beta_per_k"], float)
        sref = np.asarray(material["sips_sref"], float)
        tref = np.asarray(material["sips_tref_k"], float)
        shape = a.shape
        if shape != (len(p), 1) or any(np.asarray(v).shape != shape for v in (b, A, B, alpha, beta, sref, tref)):
            raise ValueError("competitive_sips_source requires species x one-site parameter arrays")
        if np.any(np.real(A) <= 0) or np.any(np.real(a) <= 0) or np.any(np.real(temperature) <= 0):
            raise ValueError("competitive_sips_source requires positive capacity, affinity and temperature")
        p_safe = np.where(np.real(p) > 0.0, p, 1e-30)
        p_pa = p_safe * BAR
        ninf = a * np.exp(-b / (R * temperature))
        k = A * np.exp(-B / (R * temperature))
        s = alpha * np.arctan(beta * (temperature - tref)) + sref
        if np.any(np.real(s) <= 0):
            raise ValueError("competitive_sips_source produced a non-positive Sips exponent")
        u = (k[:, 0] * p_pa) ** s[:, 0]
        den = 1.0 + np.sum(u)
        q = (ninf[:, 0] * u / den)[:, None]
        # Derivatives are with respect to engine partial pressure in bar.
        du_dp = s[:, 0] * u / p_safe
        ds_dt = alpha[:, 0] * beta[:, 0] / (1.0 + (beta[:, 0] * (temperature - tref[:, 0])) ** 2)
        dlnu_dt = ds_dt * np.log(k[:, 0] * p_pa) + s[:, 0] * B[:, 0] / (R * temperature**2)
        du_dt = u * dlnu_dt
        dninf_dt = ninf[:, 0] * b[:, 0] / (R * temperature**2)
        jac = np.empty((1, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        jac[0] = -(ninf[:, 0] * u / den**2)[:, None] * du_dp[None, :]
        jac[0][np.diag_indices(len(p))] += ninf[:, 0] * du_dp / den
        qt = (dninf_dt * u / den + ninf[:, 0] * (du_dt / den - u * np.sum(du_dt) / den**2))[:, None]
        return q, jac, qt
    if material.get("model") == "ko_extended_langmuir_source":
        # Source-labelled extended Langmuir form used by Ko (2021), Table 2.
        # The source stores pressure affinity in 1/bar and temperature
        # dependence as b=b0*exp(DE/T), with qsat=qsa+qsb/T.  This branch is
        # deliberately separate from the Ahn/Wang dual-site correlation: the
        # parameter names, units, and source material identity must not mix.
        qsa = np.asarray(material["ko_qsa_mol_kg"], float)
        qsb = np.asarray(material["ko_qsb_mol_k_kg"], float)
        b0 = np.asarray(material["ko_b0_per_bar"], float)
        de = np.asarray(material["ko_de_k"], float)
        shape = qsa.shape
        if shape != (len(p), 1) or any(np.asarray(v).shape != shape for v in (qsb, b0, de)):
            raise ValueError("ko_extended_langmuir_source requires species x one-site arrays")
        t = temperature
        qsat = qsa + qsb / t
        if np.any(np.real(qsat) <= 0.) or np.any(np.real(t) <= 0.):
            raise ValueError("Ko source equilibrium has non-positive capacity or temperature")
        p_safe = np.where(np.real(p) > 0., p, 1e-30)
        b = b0 * np.exp(de / t)
        u = b[:, 0] * p_safe
        den = 1. + np.sum(u)
        q = (qsat[:, 0] * u / den)[:, None]
        du_dp = b[:, 0]
        du_dt = u * (-de[:, 0] / t**2)
        dqsat_dt = -qsb[:, 0] / t**2
        jac = np.empty((1, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        jac[0] = -(qsat[:, 0] * u / den**2)[:, None] * du_dp[None, :]
        jac[0][np.diag_indices(len(p))] += qsat[:, 0] * du_dp / den
        qt = (dqsat_dt * u / den + qsat[:, 0] * (du_dt / den - u * np.sum(du_dt) / den**2))[:, None]
        return q, jac, qt
    if material.get("model") == "dantas_toth_source":
        # Source-local Toth branch for Dantas et al. CO2/N2--13X.  The
        # article uses independent pure-component Toth fits for the binary
        # mixture, so each component receives its own partial-pressure
        # loading and there is no Langmuir-site competition here.  K0 is in
        # bar^-1 and the stored heat is the positive magnitude of -Delta H
        # in J/mol.  This branch is deliberately not a replacement for the
        # dual-site Langmuir or Sips branches used by other sources.
        qm = np.asarray(material["dantas_qm_mol_kg"])
        exponent = np.asarray(material["dantas_n"])
        k0 = np.asarray(material["dantas_k0_per_bar"])
        heat = np.asarray(material["dantas_heat_j_mol"])
        shape = (len(p), 1)
        if any(np.asarray(value).shape != shape for value in (qm, exponent, k0, heat)):
            raise ValueError("dantas_toth_source requires species x one-site arrays")
        t = temperature
        if np.any(np.real(qm) < 0.0) or np.any(np.real(exponent) <= 0.0) or np.any(np.real(t) <= 0.0):
            raise ValueError("Dantas Toth parameters require non-negative capacity, positive exponent and temperature")
        p_safe = np.where(np.real(p) > 0.0, p, 1e-30)
        # Table 6 reports -DeltaH as a positive magnitude; Eq. (18) is
        # K=K0*exp(-DeltaH/(R*T)), hence the positive exponent below.
        k = k0 * np.exp(heat / (R * t))
        u = k * p_safe[:, None]
        un = u ** exponent
        base = 1.0 + un
        q = qm * u / base ** (1.0 / exponent)
        dq_dp = qm * k * base ** (-1.0 / exponent - 1.0)
        du_dt = -u * heat / (R * t**2)
        dq_dt = qm * base ** (-1.0 / exponent - 1.0) * du_dt
        # k=0 is used for the inert carrier in the CO2/He runs; the direct
        # derivative above already evaluates to an exact zero there.
        jac = np.zeros((1, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        jac[0][np.diag_indices(len(p))] = dq_dp[:, 0]
        return q, jac, dq_dt
    if material.get("model") == "tian_toth_source":
        # Tian 2018 source-local independent pure-component Toth fits.  This
        # branch deliberately has no competitive denominator: the article
        # reports pure CO2/CH4 fits and does not publish a binary equilibrium
        # model.  It is therefore diagnostic-only and must not be treated as
        # the source's binary VPSA model.  The source exponent is
        # n(T)=t0+c/T; retaining dn/dT matters for a non-isothermal Jacobian.
        qm = np.asarray(material["tian_qm_mol_kg"])
        t0 = np.asarray(material["tian_t0"])
        c_k = np.asarray(material["tian_c_k"])
        k0 = np.asarray(material["tian_k0_per_bar"])
        heat = np.asarray(material["tian_heat_j_mol"])
        shape = (len(p), 1)
        if any(np.asarray(value).shape != shape for value in (qm, t0, c_k, k0, heat)):
            raise ValueError("tian_toth_source requires species x one-site arrays")
        t = temperature
        if np.any(np.real(qm) < 0.0) or np.any(np.real(t) <= 0.0):
            raise ValueError("Tian Toth parameters require non-negative capacity and positive temperature")
        exponent = t0 + c_k / t
        if np.any(np.real(exponent) <= 0.0):
            raise ValueError("Tian Toth temperature-dependent exponent must be positive")
        p_safe = np.where(np.real(p) > 0.0, p, 1e-30)
        k = k0 * np.exp(heat / (R * t))
        u = k * p_safe[:, None]
        u_safe = np.where(np.real(u) > 0.0, u, 1e-300)
        un = u_safe ** exponent
        base = 1.0 + un
        q = qm * u / base ** (1.0 / exponent)
        # Exact pressure derivative after cancellation of the Toth numerator.
        dq_dp = qm * k * base ** (-1.0 / exponent - 1.0)
        # Exact temperature derivative includes both K(T) and n(T).
        dlnu_dt = -heat / (R * t**2)
        dn_dt = -c_k / t**2
        dlnbase_dt = (un / base) * (dn_dt * np.log(u_safe) + exponent * dlnu_dt)
        dlnq_dt = dlnu_dt + (dn_dt / exponent**2) * np.log(base) - dlnbase_dt / exponent
        dq_dt = q * dlnq_dt
        jac = np.zeros((1, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        jac[0][np.diag_indices(len(p))] = dq_dp[:, 0]
        return q, jac, dq_dt
    if material.get("model") == "tamof1_dual_site_langmuir_freundlich_303k":
        # TAMOF-1 source-local dual-site Langmuir--Freundlich equilibrium.
        # Supplementary Table 16 reports qm in mmol/g (numerically equal to
        # mol/kg), K in bar^(-n), and n as a dimensionless heterogeneity
        # exponent at 303 K. This branch is intentionally fixed-temperature:
        # the public source does not provide a complete CO2/CH4/N2
        # temperature-dependent parameter bundle or an LDF/thermal bundle.
        # It is therefore usable for an isolated equilibrium audit only, not
        # as an unqualified non-isothermal bed material.
        qmax = np.asarray(material["tamof_qm_mol_kg"], dtype=float)
        affinity = np.asarray(material["tamof_k_per_bar_power"], dtype=float)
        exponent = np.asarray(material["tamof_n"], dtype=float)
        if qmax.shape != (len(p), 2) or affinity.shape != qmax.shape or exponent.shape != qmax.shape:
            raise ValueError("TAMOF-1 dual-site LF requires species x two-site arrays")
        reference = float(material.get("tamof_reference_temperature_k", 303.0))
        tolerance = float(material.get("tamof_temperature_tolerance_k", 1.0e-9))
        if np.iscomplexobj(temperature) and abs(float(np.imag(temperature))) > 1.0e-12:
            raise TypeError("TAMOF-1 303 K equilibrium is not temperature-complex-step enabled")
        if abs(float(np.real(temperature)) - reference) > tolerance:
            raise ValueError(
                f"TAMOF-1 303 K equilibrium is outside its published temperature domain: "
                f"T={float(np.real(temperature)):.8g} K, reference={reference:.8g} K"
            )
        if np.any(qmax < 0.0) or np.any(affinity <= 0.0) or np.any(exponent <= 0.0):
            raise ValueError("TAMOF-1 dual-site LF parameters must be non-negative/positive")
        p_safe = np.where(np.real(p) > 0.0, p, 1.0e-30)
        u = (affinity * p_safe[:, None]) ** exponent
        den = 1.0 + np.sum(u, axis=0)
        q = qmax * u / den
        jac = np.empty((2, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        qt = np.zeros_like(q)
        for site in range(2):
            du_dp = exponent[:, site] * u[:, site] / p_safe
            jac[site] = -(qmax[:, site] * u[:, site] / den[site] ** 2)[:, None] * du_dp[None, :]
            jac[site][np.diag_indices(len(p))] += qmax[:, site] * du_dp / den[site]
        return q, jac, qt
    if material.get("model") == "chung_2025_edsl_source":
        # Shin--Yoon--Chung (2025) source-local extended dual-site
        # Langmuir equation.  The published affinity prefactor is applied on
        # a Pa pressure basis and Delta-U is stored with its reported negative
        # sign: B=b0*exp(-DeltaU/(R*T)).  Keeping these names and units
        # separate prevents accidental reuse in the historical bar-based
        # Wang/Ahn material branches.
        qsat = np.asarray(material["chung_qsat_mol_kg"], float)
        b0 = np.asarray(material["chung_b0_per_pa"], float)
        delta_u = np.asarray(material["chung_delta_u_j_mol"], float)
        shape = (len(p), 2)
        if any(array.shape != shape for array in (qsat, b0, delta_u)):
            raise ValueError("Chung EDSL requires species x two-site arrays")
        if (np.any(np.real(qsat) <= 0.0) or np.any(np.real(b0) <= 0.0)
                or np.any(np.real(delta_u) >= 0.0)
                or np.any(np.real(temperature) <= 0.0)):
            raise ValueError("Chung EDSL requires positive qsat/b0, negative Delta-U and positive temperature")
        p_safe_bar = np.where(np.real(p) > 0.0, p, 1e-30)
        p_pa = p_safe_bar * BAR
        affinity = b0 * np.exp(-delta_u / (R * temperature))
        u = affinity * p_pa[:, None]
        den = 1.0 + np.sum(u, axis=0)
        q = qsat * u / den
        du_dp_bar = affinity * BAR
        du_dt = u * delta_u / (R * temperature**2)
        dden_dt = np.sum(du_dt, axis=0)
        qt = qsat * (du_dt / den - u * dden_dt / den**2)
        jac = np.empty((2, len(p), len(p)), dtype=np.result_type(p, temperature, float))
        for site in range(2):
            diagonal = qsat[:, site] * du_dp_bar[:, site] / den[site]
            left = qsat[:, site] * u[:, site] / den[site]**2
            jac[site] = np.diag(diagonal) - left[:, None] * du_dp_bar[:, site][None, :]
        return q, jac, qt
    a1 = np.asarray(material["a1"], float)
    a2 = np.asarray(material.get("a2", np.zeros_like(a1)), float)
    energy = np.asarray(material.get("affinity_energy_j_mol", material["binding_energy_j_mol"]), float)
    b = np.asarray(material["b0_per_bar"], float) * np.exp(energy / (R * temperature))
    cap = a1 - a2 * temperature
    check_capacity(material,temperature,cap)
    den = 1.0 + np.sum(b * p[:, None], axis=0)
    q = cap * b * p[:, None] / den
    db = -b * energy / (R * temperature**2)
    dden = np.sum(db * p[:, None], axis=0)
    qt = -a2 * b * p[:, None] / den + cap * p[:, None] * (db / den - b * dden / den**2)
    jac = np.empty((a1.shape[1], len(p), len(p)),dtype=np.result_type(p,temperature,float))
    for site in range(a1.shape[1]):
        jac[site] = np.diag(cap[:, site] * b[:, site] / den[site]) - (
            (cap[:, site] * b[:, site] * p / den[site]**2)[:, None] * b[:, site][None, :]
        )
    return q, jac, qt


def valve_molar_flow(p_up_bar, p_down_bar, temperature, y, mw, cp, area_cd_m2, regularization_bar=1e-6):
    """One-way ideal-gas compressible orifice, including critical choking.

    area_cd is discharge coefficient times area, not a fitted mol/s/bar gain.
    A closed or reverse pressure gradient gives exactly zero flow.
    """
    if np.real(p_up_bar) <= np.real(p_down_bar) or area_cd_m2 <= 0:
        return 0.0
    mixture_mw = np.dot(y, mw)
    cp_mix = np.dot(y, cp)
    gamma = cp_mix / (cp_mix - R)
    ratio = _positive(p_down_bar / p_up_bar)
    critical = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    if np.real(ratio) <= np.real(critical):
        psi = np.sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
    else:
        # Algebraically identical to subtracting the two powers of p_down/
        # p_up, but accurate when they are almost equal. No valve-law fitting.
        log_ratio=np.log1p(-(p_up_bar-p_down_bar)/p_up_bar)
        psi = np.sqrt(2*gamma/(gamma-1)*np.exp(2/gamma*log_ratio)*
                      (-np.expm1((gamma-1)/gamma*log_ratio)))
    flow=area_cd_m2 * p_up_bar * BAR * psi / np.sqrt(mixture_mw * R * temperature)
    # Physical valves have viscous behavior near zero driving pressure. The
    # inviscid sqrt law has an infinite derivative and chatters after PE.
    delta=p_up_bar-p_down_bar
    if np.real(delta)<regularization_bar:
        flow*=np.sqrt(delta/regularization_bar)
    return flow


class Bed:
    """A bed with local gas/solid thermal equilibrium and an independent wall."""

    def __init__(self, spec, species):
        self.spec = spec
        self.ns = len(species)
        self.n = int(spec["n_cells"])
        self.nv = 2 * self.ns + 2
        self.size = self.n * self.nv
        self.cp = np.array([s["cp_j_mol_k"] for s in species])
        self.cv = self.cp - R
        self.mw = np.array([s["mw_kg_mol"] for s in species])
        self.equation_of_state = spec.get("equation_of_state", "ideal_gas")
        self.pr_components = None
        if self.equation_of_state == "peng_robinson":
            required = ("critical_temperature_k", "critical_pressure_pa", "acentric_factor")
            missing = [f"{entry.get('name', index)}:{field}" for index, entry in enumerate(species)
                       for field in required if field not in entry]
            if missing:
                raise ValueError(f"Peng-Robinson species properties are missing: {missing}")
            interaction = spec.get("pr_binary_interaction")
            self.pr_components = PengRobinsonComponents(
                names=tuple(str(entry["name"]) for entry in species),
                critical_temperature_k=tuple(float(entry["critical_temperature_k"]) for entry in species),
                critical_pressure_pa=tuple(float(entry["critical_pressure_pa"]) for entry in species),
                acentric_factor=tuple(float(entry["acentric_factor"]) for entry in species),
                binary_interaction=None if interaction is None else tuple(tuple(float(value) for value in row) for row in interaction),
            )
            # Trigger all dimension, symmetry and positivity checks at setup.
            self.pr_components.arrays()
        elif self.equation_of_state != "ideal_gas":
            raise ValueError(f"Unsupported equation of state: {self.equation_of_state}")
        self.pr_departure_energy = bool(spec.get("pr_departure_energy", False))
        if self.pr_departure_energy and self.equation_of_state != "peng_robinson":
            raise ValueError("Peng-Robinson departure energy requires equation_of_state='peng_robinson'")
        self.area = math.pi * spec["diameter_m"]**2 / 4
        self.materials = spec["materials"]
        ends=np.array([m["end_fraction"] for m in self.materials],float)
        fractions=np.diff(np.r_[0.,ends])
        if len(fractions)>self.n or np.any(fractions<=0) or abs(ends[-1]-1)>1e-10:
            raise ValueError("Layer fractions must increase to one with at least one cell per layer")
        # Allocate at least one cell per layer and preserve its exact length.
        # Midpoint material selection on a uniform grid changes the adsorbent
        # ratio whenever a layer boundary is not an integer cell boundary.
        counts=np.ones(len(fractions),dtype=int)
        for _ in range(self.n-len(fractions)):
            counts[np.argmax(fractions/(counts+1))]+=1
        self.layer_indices=[material for material,count in zip(self.materials,counts) for _ in range(count)]
        self.dz=np.concatenate([np.full(count,fraction*spec["length_m"]/count) for fraction,count in zip(fractions,counts)])
        self.cell_centers_m=np.cumsum(self.dz)-self.dz/2
        self.face_spacing=(self.dz[:-1]+self.dz[1:])/2
        self.cell_volume = self.area * self.dz
        # Optional layer-specific hydraulic/thermal properties are used by
        # source-labelled layered-bed screens.  Existing cases keep their
        # historical bed-global defaults when a material does not override a
        # field, so this extension does not silently change S-015.
        self.void_fraction = np.asarray([
            float(material.get("void_fraction", spec["void_fraction"]))
            for material in self.layer_indices
        ])
        self.pellet_void_fraction = np.asarray([
            float(material.get("pellet_void_fraction", spec.get("pellet_void_fraction", 0.0)))
            for material in self.layer_indices
        ])
        self.particle_diameter = np.asarray([
            float(material.get("particle_diameter_m", spec["particle_diameter_m"]))
            for material in self.layer_indices
        ])
        self.bulk_density = np.asarray([
            float(material.get("bulk_density_kg_m3", spec["bulk_density_kg_m3"]))
            for material in self.layer_indices
        ])
        self.solid_cp = np.asarray([
            float(material.get("solid_cp_j_kg_k", spec["solid_cp_j_kg_k"]))
            for material in self.layer_indices
        ])
        # Source-specific thermal closures may define the solid/adsorbate
        # storage basis separately from the mass-balance basis.  The default
        # factor is one, preserving all historical H2/Wang behavior.  A
        # Büchner--Pini case can explicitly bind the printed rho_b/epsilon
        # energy basis without corrupting gas inventory or pressure drop.
        self.energy_capacity_factor = np.asarray([
            float(material.get("energy_capacity_factor", spec.get("energy_capacity_factor", 1.0)))
            for material in self.layer_indices
        ])
        self.conduction_factor = np.asarray([
            float(material.get("conduction_factor", spec.get("conduction_factor", 1.0)))
            for material in self.layer_indices
        ])
        self.bed_k = np.asarray([
            float(material.get("bed_k_w_m_k", spec.get("bed_k_w_m_k", 0.24)))
            for material in self.layer_indices
        ])
        self.vgas = self.void_fraction * self.cell_volume
        self.vparticle = (1 - self.void_fraction) * self.cell_volume
        self.mass = self.bulk_density * self.cell_volume
        self.cs = self.mass * self.solid_cp
        self.thermal_mass = self.mass * self.energy_capacity_factor
        self.thermal_cs = self.cs * self.energy_capacity_factor
        self.alpha = self.vparticle * self.pellet_void_fraction * BAR / R
        self.wall_area = math.pi * spec["diameter_m"] * self.dz
        radius = spec["diameter_m"] / 2
        outer_radius = radius + spec.get("wall_thickness_m", 0.002)
        self.wall_cross_area = math.pi * (outer_radius**2 - radius**2)
        self.outer_area = 2 * math.pi * outer_radius * self.dz
        # Source papers do not always report the same wall-capacity basis.
        # Accept an explicitly volumetric heat capacity, but never guess that
        # a printed J/kg/K value really meant J/m3/K.  Reject a mixed basis
        # rather than silently multiplying the same material capacity twice.
        volumetric_wall_cp = spec.get("wall_volumetric_cp_j_m3_k")
        if volumetric_wall_cp is not None:
            if "wall_density_kg_m3" in spec or "wall_cp_j_kg_k" in spec:
                raise ValueError("Specify wall heat capacity as either volumetric or density times specific heat, not both")
            wall_capacity_per_volume = float(volumetric_wall_cp)
            self.wall_capacity_basis = "explicit_volumetric_j_m3_k"
        else:
            wall_capacity_per_volume = float(spec.get("wall_density_kg_m3", 7850)) * float(spec.get("wall_cp_j_kg_k", 460))
            self.wall_capacity_basis = "density_times_specific_heat"
        if not math.isfinite(wall_capacity_per_volume) or wall_capacity_per_volume <= 0:
            raise ValueError("Wall volumetric heat capacity must be finite and positive")
        self.cwall = self.wall_cross_area * self.dz * wall_capacity_per_volume
        self.wall_model = spec.get("wall_model", "dynamic_lumped")
        if self.wall_model not in {"dynamic_lumped", "fixed_ambient"}:
            raise ValueError("wall_model must be 'dynamic_lumped' or 'fixed_ambient'")
        self.material_batches=[]
        offset=0
        for material,count in zip(self.materials,counts):
            self.material_batches.append((slice(offset,offset+count),material))
            offset+=count
        self.mode = spec.get("kinetics", "loading_ldf")
        self.energy_mode = spec.get("energy_mode", "internal_energy")
        if self.pr_departure_energy and self.mode != "loading_ldf":
            raise ValueError("PR departure-energy closure currently requires loading_ldf gas states")
        if self.pr_departure_energy and self.energy_mode != "internal_energy":
            raise ValueError("PR departure-energy closure requires internal_energy mode")
        # Optional axial mass dispersion term.  The historical H2/Wang
        # cases default to zero, preserving their original upwind transport;
        # source-labelled breakthrough adapters may bind a measured or
        # independently estimated coefficient in m2/s.
        self.axial_dispersion_m2_s = float(spec.get("axial_dispersion_m2_s", 0.0))
        self.axial_dispersion_model = spec.get("axial_dispersion_model")
        if not math.isfinite(self.axial_dispersion_m2_s) or self.axial_dispersion_m2_s < 0:
            raise ValueError("axial_dispersion_m2_s must be finite and nonnegative")
        if self.axial_dispersion_model not in {None, "guderian_2021_eq3"}:
            raise ValueError("Unsupported axial_dispersion_model")
        if (self.axial_dispersion_model == "guderian_2021_eq3"
                and (tuple(entry["name"] for entry in species) != ("O2", "N2")
                     or self.ns != 2)):
            raise ValueError("Guderian Eq. 3 dispersion requires O2/N2 source order")
        if self.axial_dispersion_model is not None and self.axial_dispersion_m2_s != 0.0:
            raise ValueError("Use either a dynamic axial-dispersion model or a constant, not both")
        self.dispersion_particle_radius_m = spec.get("dispersion_particle_radius_m")
        if self.axial_dispersion_model == "guderian_2021_eq3":
            if (self.dispersion_particle_radius_m is None
                    or not math.isfinite(float(self.dispersion_particle_radius_m))
                    or float(self.dispersion_particle_radius_m) <= 0.0):
                raise ValueError("Guderian Eq. 3 requires explicit dispersion_particle_radius_m")
            self.dispersion_particle_radius_m = float(self.dispersion_particle_radius_m)
        if self.n < 2 or self.mode not in {"loading_ldf","partial_pressure_ldf"}:
            raise ValueError("Need >=2 finite volumes and a supported kinetic law")
        if self.energy_mode not in {"internal_energy","source_sensible"}:
            raise ValueError("Unknown thermal closure")
        if (not 0 < spec["void_fraction"] < 1 or np.any(self.void_fraction <= 0)
                or np.any(self.void_fraction >= 1) or np.any(self.particle_diameter <= 0)
                or np.any(self.mass <= 0) or np.any(self.cs <= 0)
                or np.any(self.energy_capacity_factor <= 0)
                or np.any(self.conduction_factor <= 0)
                or np.any(self.cwall <= 0)):
            raise ValueError("Invalid bed porosity or heat capacity")
        for material in self.materials:
            model = material.get("model", "dual_site_langmuir")
            if model == "competitive_sips_source":
                shape=np.asarray(material["sips_a_mol_kg"]).shape
                fields=("sips_b_j_mol", "sips_A_per_pa", "sips_B_j_mol", "sips_alpha",
                        "sips_beta_per_k", "sips_sref", "sips_tref_k", "binding_energy_j_mol")
            elif model == "extended_langmuir_3":
                shape=np.asarray(material["ip1_kmol_kg"]).shape
                fields=("ip2_kmol_kg_k", "ip3_per_bar", "ip4_k", "binding_energy_j_mol")
            elif model == "extended_langmuir_freundlich":
                shape=np.asarray(material["lf_k1_mol_kg"]).shape
                fields=("lf_k2_mol_kg_k","lf_k3_per_atm_power","lf_k4_k",
                        "lf_k5","lf_k6_k","binding_energy_j_mol")
            elif model == "dantas_toth_source":
                # Dantas uses one independent Toth fit per component.  Keep
                # this source-local contract separate from competitive
                # Langmuir/Sips tables used by the other gas systems.
                shape=np.asarray(material["dantas_qm_mol_kg"]).shape
                fields=("dantas_n", "dantas_k0_per_bar", "dantas_heat_j_mol",
                        "binding_energy_j_mol")
            elif model == "tian_toth_source":
                # Tian's pure-component Toth namespace is intentionally
                # separate from Dantas and every competitive source branch.
                shape=np.asarray(material["tian_qm_mol_kg"]).shape
                fields=("tian_t0", "tian_c_k", "tian_k0_per_bar", "tian_heat_j_mol",
                        "binding_energy_j_mol")
            elif model == "tamof1_dual_site_langmuir_freundlich_303k":
                # TAMOF-1 Table-16 equilibrium is a source-local dual-site
                # Langmuir--Freundlich namespace.  It is deliberately kept
                # separate from the historical dual-site Langmuir arrays:
                # the source publishes no LDF or complete thermal bundle.
                # Zero binding energies are accepted only for the explicitly
                # isothermal diagnostic case factory; they do not authorize a
                # non-isothermal or cyclic release.
                shape=np.asarray(material["tamof_qm_mol_kg"]).shape
                fields=("tamof_k_per_bar_power", "tamof_n")
            elif model == "chung_2025_edsl_source":
                # Source Table S1 uses Pa-based affinity prefactors and
                # signed internal-energy changes.  Do not alias these arrays
                # to the generic bar-based dual-site material contract.
                shape=np.asarray(material["chung_qsat_mol_kg"]).shape
                fields=("chung_b0_per_pa", "chung_delta_u_j_mol",
                        "binding_energy_j_mol")
            elif model == "streb_sips_wilson_rast_candidate":
                if self.mode != "loading_ldf":
                    raise ValueError("Streb Wilson/RAST candidate requires loading_ldf kinetics")
                if tuple(material.get("species_order", ())) != tuple(s["name"] for s in species):
                    raise ValueError("Streb Wilson/RAST species order must match the case")
                if set(material["species_order"]) != {"H2", "CO2", "CH4"} or self.ns != 3:
                    raise ValueError("Streb Wilson/RAST candidate requires H2, CO2 and CH4 exactly")
                from .streb_2021_wilson import _validated_matrix
                lam = _validated_matrix(material.get("wilson_lambda_matrix"), self.ns)
                h2 = material["species_order"].index("H2")
                if not (np.allclose(lam[h2, :], 1.0, rtol=0, atol=1e-12)
                        and np.allclose(lam[:, h2], 1.0, rtol=0, atol=1e-12)):
                    raise ValueError("Streb Wilson/RAST candidate requires unity H2 interaction pairs")
                if not isinstance(material.get("parameter_provenance"), str) or not material["parameter_provenance"].strip():
                    raise ValueError("Streb Wilson/RAST explicit matrix provenance is required")
                shape=(self.ns, 1)
                fields=("binding_energy_j_mol",)
            elif model == "guderian_2021_sips_iast":
                if self.mode != "loading_ldf":
                    raise ValueError("Guderian Sips/IAST/QDF requires loading-state kinetics")
                if tuple(material.get("species_order", ())) != tuple(s["name"] for s in species):
                    raise ValueError("Guderian source species order must match the case")
                if tuple(material.get("species_order", ())) != ("O2", "N2") or self.ns != 2:
                    raise ValueError("Guderian source model requires O2 and N2 exactly")
                if material.get("kinetic_model") != "guderian_qdf_darken":
                    raise ValueError("Guderian source model requires its published Darken-QDF kinetics")
                provenance = material.get("parameter_provenance")
                if not isinstance(provenance, str) or not provenance.strip():
                    raise ValueError("Guderian source parameter provenance is required")
                shape=(self.ns, 1)
                fields=()
            else:
                shape=np.asarray(material["a1"]).shape
                fields=("b0_per_bar","binding_energy_j_mol")
            if len(shape)!=2 or shape[0]!=self.ns:
                raise ValueError("Each material table must have species x adsorption-site axes")
            if any(np.asarray(material[k]).shape!=shape for k in fields):
                raise ValueError("Inconsistent adsorption parameter dimensions")
            binding_shape = shape
            if model == "tamof1_dual_site_langmuir_freundlich_303k" and self.mode == "loading_ldf":
                # Loading-LDF stores one total loading per species while the
                # source equilibrium remains a two-site table.  The zero
                # binding-energy diagnostic is consequently species x one
                # site, unlike the two-site equilibrium coefficients.
                binding_shape = (self.ns, 1)
            if (model != "guderian_2021_sips_iast"
                    and np.asarray(material["binding_energy_j_mol"]).shape != binding_shape):
                raise ValueError("Inconsistent binding-energy dimensions")
            if model != "extended_langmuir_freundlich" and "affinity_energy_j_mol" in material and np.asarray(material["affinity_energy_j_mol"]).shape!=shape:
                raise ValueError("Inconsistent affinity temperature coefficient dimensions")
            if model != "guderian_2021_sips_iast":
                if np.asarray(material["ldf_per_s"]).shape!=(self.ns,) or min(material["ldf_per_s"])<0:
                    raise ValueError("Invalid kinetic coefficient vector")
        if self.mode == "partial_pressure_ldf" and np.any(self.alpha <= 0):
            raise ValueError("Explicit pellet gas requires a positive pellet void fraction")

    def initial(self, p, y, temperature):
        x = np.zeros((self.n, self.nv))
        y = np.asarray(y, dtype=float)
        if self.equation_of_state == "peng_robinson":
            density = pr_molar_density_from_pressure(p * BAR, temperature, y, self.pr_components)
            x[:, :self.ns] = y[None, :] * density * self.vgas[:, None]
        else:
            x[:, :self.ns] = y[None,:] * p * BAR * self.vgas[:,None] / (R * temperature)
        for j, material in enumerate(self.layer_indices):
            q, _, _ = equilibrium(p * y, temperature, material)
            x[j, self.ns:2*self.ns] = p * y if self.mode == "partial_pressure_ldf" else q.sum(axis=1)
        x[:, -2:] = temperature
        return x.ravel()

    def fields(self, state):
        x = np.asarray(state).reshape(self.n, self.nv)
        gas = _positive(x[:, :self.ns],1e-30)
        ads = _positive(x[:, self.ns:2*self.ns])
        t, tw = x[:, -2], x[:, -1]
        if np.any(np.real(t) <= 100) or np.any(np.real(t) >= 1000) or np.any(np.real(tw) <= 100):
            raise FloatingPointError("Bed temperature outside physical numerical domain")
        total = gas.sum(axis=1)
        fractions = gas / total[:, None]
        if self.equation_of_state == "peng_robinson":
            pressure = np.asarray([
                pr_pressure_from_molar_density(total[index] / self.vgas[index], t[index], fractions[index], self.pr_components)
                for index in range(self.n)
            ]) / BAR
        else:
            pressure = total * R * t / self.vgas / BAR
        return gas, ads, t, tw, pressure, fractions

    def inventory_energy(self, state):
        gas, ads, t, tw, _, _ = self.fields(state)
        inventory = gas.sum(axis=0)
        energy = np.sum(self.thermal_cs * (t - 298.15))
        if self.wall_model == "dynamic_lumped":
            energy += np.sum(self.cwall * (tw - 298.15))
        for j, material in enumerate(self.layer_indices):
            model = material.get("model", "dual_site_langmuir")
            if model == "guderian_2021_sips_iast":
                from .guderian_2021_n2 import (
                    SPECIES_ORDER,
                    adsorption_energy_magnitude_j_kg,
                )
                inventory += self.mass[j] * ads[j]
                stored = sum(
                    adsorption_energy_magnitude_j_kg(name, ads[j, index], t[j])
                    for index, name in enumerate(SPECIES_ORDER)
                )
                energy -= self.thermal_mass[j] * stored
                if self.energy_mode == "internal_energy":
                    energy += self.gas_internal_energy(gas[j], t[j], j)
                continue
            weights = np.asarray(material["binding_energy_j_mol"])
            if self.mode == "partial_pressure_ldf":
                q, _, _ = equilibrium(ads[j], t[j], material)
                pore = self.alpha[j] * ads[j] / t[j]
                inventory += pore + self.mass[j] * q.sum(axis=1)
                energy -= self.thermal_mass[j] * np.sum(weights * q)
                if self.energy_mode == "internal_energy":
                    energy += np.dot(gas[j] + pore, self.cv) * t[j]
                    if self.pr_departure_energy:
                        energy += pr_total_internal_energy_departure(
                            gas[j], self.vgas[j], t[j], self.pr_components
                        )
            else:
                if weights.shape[1] != 1:
                    raise ValueError("Loading LDF requires a single-site binding energy")
                inventory += self.mass[j] * ads[j]
                energy -= self.thermal_mass[j] * np.dot(weights[:, 0], ads[j])
                if self.energy_mode == "internal_energy":
                    energy += self.gas_internal_energy(gas[j], t[j], j)
        return inventory, float(energy)

    def gas_internal_energy(self, inventory, temperature, cell_index):
        energy = np.dot(inventory, self.cv) * temperature
        if self.pr_departure_energy:
            energy += pr_total_internal_energy_departure(
                inventory, self.vgas[cell_index], temperature, self.pr_components
            )
        return energy

    def gas_internal_energy_derivatives(self, inventory, temperature, cell_index):
        component_derivative = self.cv * temperature
        temperature_derivative = np.dot(inventory, self.cv)
        if self.pr_departure_energy:
            departure_component, departure_temperature = pr_total_internal_energy_departure_derivatives(
                np.real(inventory), self.vgas[cell_index], float(np.real(temperature)), self.pr_components
            )
            component_derivative = component_derivative + departure_component
            temperature_derivative += departure_temperature
        return component_derivative, temperature_derivative

    def enthalpy(self, y, t, pressure_bar=None):
        reference = 0.0 if self.energy_mode == "internal_energy" else self.spec["ambient_k"]
        enthalpy = np.dot(y, self.cp) * (t - reference)
        if self.pr_departure_energy:
            if pressure_bar is None:
                raise ValueError("PR departure enthalpy requires the stream pressure")
            enthalpy += pr_molar_enthalpy_departure(
                float(np.real(pressure_bar)) * BAR,
                float(np.real(t)),
                np.real(y),
                self.pr_components,
            )
        return enthalpy

    def internal_flux(self, state):
        return vectorized.internal_flux(self,state)

    def internal_flux_scalar(self, state):
        """Preserved scalar reference for independent vectorization tests."""
        gas, _, t, _, p, y = self.fields(state)
        dtype=np.result_type(state,float)
        flux = np.zeros((self.n + 1, self.ns),dtype=dtype)
        heat = np.zeros(self.n + 1,dtype=dtype)
        eps = 0.5 * (self.void_fraction[:-1] + self.void_fraction[1:])
        dpellet = 0.5 * (self.particle_diameter[:-1] + self.particle_diameter[1:])
        linear = 150 * self.spec.get("viscosity_pa_s", 1.2e-5) * (1-eps)**2 / (eps**3 * dpellet**2)
        for face in range(1, self.n):
            # ``eps``, ``dpellet`` and ``linear`` are face arrays.  The
            # scalar reference must select the current face explicitly;
            # otherwise the vectorized expressions below broadcast against
            # the component vector and corrupt the scalar flux shape.
            face_index = face - 1
            eps_face = eps[face_index]
            dpellet_face = dpellet[face_index]
            linear_face = linear[face_index]
            delta = (p[face-1] - p[face]) * BAR
            pf = (p[face-1] + p[face]) * BAR / 2
            tf = (t[face-1] + t[face]) / 2
            yf = (y[face-1] + y[face]) / 2
            if self.equation_of_state == "peng_robinson":
                concentration = 0.5 * (gas[face-1].sum() / self.vgas[face-1]
                                       + gas[face].sum() / self.vgas[face])
                rho = concentration * np.dot(yf, self.mw)
            else:
                rho = pf * np.dot(yf, self.mw) / (R * tf)
            quadratic = 1.75 * rho * (1-eps_face) / (eps_face**3 * dpellet_face)
            sign=1 if np.real(delta)>=0 else -1
            gradient = sign*delta / self.face_spacing[face-1]
            # This form avoids subtractive cancellation in the Darcy limit.
            velocity = 2 * gradient / (linear_face + np.sqrt(linear_face**2 + 4*quadratic*gradient))
            if self.equation_of_state == "peng_robinson":
                flow = sign * velocity * self.area * concentration
            else:
                flow = sign*velocity * self.area * pf / (R * tf)
            donor = face-1 if sign > 0 else face
            flux[face] = flow * y[donor]
            heat[face] = flow * self.enthalpy(y[donor], t[donor], p[donor])
            if self.axial_dispersion_model == "guderian_2021_eq3":
                from .guderian_2021_n2 import axial_dispersion_o2_n2_m2_s
                dispersion = axial_dispersion_o2_n2_m2_s(
                    velocity, pf / BAR, tf, eps_face, self.dispersion_particle_radius_m
                )
                concentration = 0.5 * (gas[face-1].sum() / self.vgas[face-1]
                                       + gas[face].sum() / self.vgas[face])
                diffusive = (-dispersion * self.area * concentration
                              * (y[face] - y[face-1]) / self.face_spacing[face-1])
                flux[face] += diffusive
                yf = 0.5 * (y[face-1] + y[face])
                tf = 0.5 * (t[face-1] + t[face])
                heat[face] += diffusive.sum() * self.enthalpy(yf, tf, pf / BAR)
            elif self.axial_dispersion_m2_s > 0:
                concentration = 0.5 * (gas[face-1].sum() / self.vgas[face-1]
                                       + gas[face].sum() / self.vgas[face])
                diffusive = (-self.axial_dispersion_m2_s * self.area * concentration
                              * (y[face] - y[face-1]) / self.face_spacing[face-1])
                flux[face] += diffusive
                yf = 0.5 * (y[face-1] + y[face])
                tf = 0.5 * (t[face-1] + t[face])
                heat[face] += diffusive.sum() * self.enthalpy(yf, tf, pf / BAR)
        return flux, heat

    def rates(self, state, flux, enthalpy_flux, h_bed_wall_override=None):
        return vectorized.rates(self,state,flux,enthalpy_flux,
                                h_bed_wall_override=h_bed_wall_override)

    def rates_scalar(self, state, flux, enthalpy_flux, h_bed_wall_override=None):
        """Preserved scalar reference; physical equations must stay identical."""
        gas, ads, t, tw, p, y = self.fields(state)
        flow_rate = flux[:-1] - flux[1:]
        heat = enthalpy_flux[:-1] - enthalpy_flux[1:]
        dtype=np.result_type(state,float)
        qbed = np.zeros(self.n,dtype=dtype)
        qwall = np.zeros(self.n,dtype=dtype)
        gb = (0.5 * (self.bed_k[:-1] * self.conduction_factor[:-1]
                     + self.bed_k[1:] * self.conduction_factor[1:])) * self.area / self.face_spacing
        gw = self.spec.get("wall_k_w_m_k", 10.9) * self.wall_cross_area / self.face_spacing
        qb = gb * np.diff(t)
        qw = gw * np.diff(tw)
        qbed[:-1] += qb; qbed[1:] -= qb
        qwall[:-1] += qw; qwall[1:] -= qw
        h_bed_wall = (self.spec.get("h_bed_wall_w_m2_k", 0) if h_bed_wall_override is None
                      else h_bed_wall_override)
        if self.wall_model == "fixed_ambient":
            wall_exchange = h_bed_wall * self.wall_area * (self.spec["ambient_k"]-t)
            ambient = wall_exchange.copy()
        else:
            wall_exchange = h_bed_wall * self.wall_area * (tw-t)
            ambient = self.spec.get("h_ambient_w_m2_k", 0) * self.outer_area * (self.spec["ambient_k"]-tw)
        dx = np.zeros((self.n, self.nv),dtype=dtype)
        if self.wall_model == "fixed_ambient":
            dx[:, -1] = 0.0
        else:
            dx[:, -1] = (-wall_exchange + qwall + ambient) / self.cwall
        for j, material in enumerate(self.layer_indices):
            model = material.get("model", "dual_site_langmuir")
            if self.mode == "loading_ldf":
                if model == "guderian_2021_sips_iast":
                    q, _, _ = equilibrium(p[j]*y[j], t[j], material)
                    from .guderian_2021_n2 import (
                        SPECIES_ORDER,
                        adsorption_energy_temperature_derivative_j_kg_k,
                        heat_of_adsorption_magnitude_j_mol,
                        qdf_rate_mol_kg_s,
                    )
                    da = qdf_rate_mol_kg_s(ads[j], q.sum(axis=1), t[j])
                    heat_weights = np.asarray([
                        heat_of_adsorption_magnitude_j_mol(name, ads[j, index], t[j])
                        for index, name in enumerate(SPECIES_ORDER)
                    ])
                    energy_temperature_derivative = sum(
                        adsorption_energy_temperature_derivative_j_kg_k(
                            name, ads[j, index], t[j]
                        )
                        for index, name in enumerate(SPECIES_ORDER)
                    )
                else:
                    ldf = np.asarray(material["ldf_per_s"])
                    heat_weights = np.asarray(material["binding_energy_j_mol"])[:, 0]
                    energy_temperature_derivative = 0.0
                    # See the vectorized implementation: an all-zero LDF is
                    # an exact inert-material contract, so evaluating an
                    # unused equilibrium expression is both unnecessary and
                    # numerically unsafe during implicit-solver trial states.
                    if np.all(ldf == 0.0):
                        da = np.zeros_like(ads[j])
                    else:
                        q, _, _ = equilibrium(p[j]*y[j], t[j], material)
                        da = ldf * (q.sum(axis=1)-ads[j])
                dn = flow_rate[j] - self.mass[j] * da
                rhs = heat[j] + qbed[j] + wall_exchange[j] + self.thermal_mass[j] * np.dot(heat_weights, da)
                capacity = self.thermal_cs[j] - self.thermal_mass[j] * energy_temperature_derivative
                if self.energy_mode == "internal_energy":
                    gas_component_energy, gas_temperature_capacity = self.gas_internal_energy_derivatives(
                        gas[j], t[j], j
                    )
                    rhs -= np.dot(dn, gas_component_energy)
                    capacity += gas_temperature_capacity
                if np.real(capacity) <= 0.0:
                    raise FloatingPointError("Guderian/source thermal capacity became non-positive")
                dt = rhs / capacity
            else:
                weights = np.asarray(material["binding_energy_j_mol"])
                ldf = np.asarray(material["ldf_per_s"])
                q, jac_site, qt_site = equilibrium(ads[j], t[j], material)
                jac = jac_site.sum(axis=0)
                qt = qt_site.sum(axis=1)
                exchange = self.vparticle[j] * ldf * (p[j]*y[j]-ads[j]) * BAR / (R*t[j])
                dn = flow_rate[j] - exchange
                matrix = np.zeros((self.ns+1, self.ns+1),dtype=dtype)
                matrix[:self.ns, :self.ns] = self.alpha[j]/t[j]*np.eye(self.ns) + self.mass[j]*jac
                matrix[:self.ns, -1] = -self.alpha[j]*ads[j]/t[j]**2 + self.mass[j]*qt
                binding_p = sum(weights[:, s] @ jac_site[s] for s in range(weights.shape[1]))
                binding_t = np.sum(weights*qt_site)
                matrix[-1, :self.ns] = -self.thermal_mass[j] * binding_p
                matrix[-1, -1] = self.thermal_cs[j] - self.thermal_mass[j] * binding_t
                rhs_heat = heat[j] + qbed[j] + wall_exchange[j]
                if self.energy_mode == "internal_energy":
                    matrix[-1, :self.ns] += self.alpha[j]*self.cv
                    matrix[-1, -1] += np.dot(gas[j], self.cv)
                    rhs_heat -= np.dot(dn, self.cv)*t[j]
                solved = np.linalg.solve(matrix, np.r_[exchange, rhs_heat])
                da, dt = solved[:-1], solved[-1]
            dx[j, :self.ns] = dn
            dx[j, self.ns:2*self.ns] = da
            dx[j, -2] = dt
        return dx.ravel(), ambient.sum()


class ProductReceiver:
    """Rigid, perfectly mixed product receiver with conservative inventories.

    The receiver is an optional *physical* state, not a cumulative stream
    ledger.  It stores one molar inventory per component plus one gas
    temperature.  The same gas equation of state and sensible internal-energy
    convention as the connected bed are used.  This is required by source
    cycles in which product backflow reverses into a depressurised bed and the
    receiver simultaneously supplies purge and an external product stream.

    When the connected bed explicitly enables ``pr_departure_energy``, the
    receiver uses the same residual internal-energy and stream-enthalpy
    closure.  Otherwise historical ideal-gas sensible energy is preserved.
    """

    def __init__(self, spec, bed):
        self.spec = dict(spec)
        self.ns = bed.ns
        self.cp = bed.cp.copy()
        self.cv = bed.cv.copy()
        self.equation_of_state = bed.equation_of_state
        self.pr_components = bed.pr_components
        self.pr_departure_energy = bed.pr_departure_energy
        self.volume_m3 = float(spec["volume_m3"])
        self.thermal_mode = spec.get("thermal_mode", "adiabatic_sensible")
        if not np.isfinite(self.volume_m3) or self.volume_m3 <= 0.0:
            raise ValueError("Product receiver volume_m3 must be finite and positive")
        if self.thermal_mode != "adiabatic_sensible":
            raise ValueError("Product receiver thermal_mode must be 'adiabatic_sensible'")
        self.size = self.ns + 1

    def initial(self, pressure_bar, composition, temperature_k):
        composition = np.asarray(composition, dtype=float)
        if (composition.shape != (self.ns,) or np.any(~np.isfinite(composition))
                or np.any(composition < 0.0) or not np.isclose(composition.sum(), 1.0, atol=1e-10)):
            raise ValueError("Product receiver composition must be finite, nonnegative and sum to one")
        pressure_bar = float(pressure_bar)
        temperature_k = float(temperature_k)
        if not np.isfinite(pressure_bar) or pressure_bar <= 0.0:
            raise ValueError("Product receiver initial pressure must be finite and positive")
        if not np.isfinite(temperature_k) or not 100.0 < temperature_k < 1000.0:
            raise ValueError("Product receiver initial temperature is outside the numerical domain")
        if self.equation_of_state == "peng_robinson":
            density = pr_molar_density_from_pressure(
                pressure_bar * BAR, temperature_k, composition, self.pr_components
            )
        else:
            density = pressure_bar * BAR / (R * temperature_k)
        return np.r_[composition * density * self.volume_m3, temperature_k]

    def fields(self, state):
        state = np.asarray(state)
        if state.shape != (self.size,):
            raise ValueError("Product receiver state has incompatible dimensions")
        inventory = _positive(state[:self.ns], 1e-30)
        temperature = state[-1]
        if np.real(temperature) <= 100.0 or np.real(temperature) >= 1000.0:
            raise FloatingPointError("Product receiver temperature outside physical numerical domain")
        total = inventory.sum()
        composition = inventory / total
        density = total / self.volume_m3
        if self.equation_of_state == "peng_robinson":
            pressure = pr_pressure_from_molar_density(
                density, temperature, composition, self.pr_components
            ) / BAR
        else:
            pressure = density * R * temperature / BAR
        return inventory, temperature, pressure, composition

    def enthalpy(self, composition, temperature, pressure_bar=None):
        enthalpy = np.dot(composition, self.cp) * temperature
        if self.pr_departure_energy:
            if pressure_bar is None:
                raise ValueError("PR receiver departure enthalpy requires pressure")
            enthalpy += pr_molar_enthalpy_departure(
                float(np.real(pressure_bar)) * BAR,
                float(np.real(temperature)),
                np.real(composition),
                self.pr_components,
            )
        return enthalpy

    def inventory_energy(self, state):
        inventory, temperature, _, _ = self.fields(state)
        energy = np.dot(inventory, self.cv) * temperature
        if self.pr_departure_energy:
            energy += pr_total_internal_energy_departure(
                inventory, self.volume_m3, temperature, self.pr_components
            )
        return inventory.copy(), float(np.real(energy))

    def rates(self, state, component_rate, enthalpy_rate_w):
        inventory, temperature, _, _ = self.fields(state)
        component_rate = np.asarray(component_rate)
        if component_rate.shape != (self.ns,):
            raise ValueError("Product receiver component rate has incompatible dimensions")
        component_energy = self.cv * temperature
        capacity = np.dot(inventory, self.cv)
        if self.pr_departure_energy:
            departure_component, departure_temperature = pr_total_internal_energy_departure_derivatives(
                np.real(inventory), self.volume_m3, float(np.real(temperature)), self.pr_components
            )
            component_energy = component_energy + departure_component
            capacity = capacity + departure_temperature
        if np.real(capacity) <= 0.0:
            raise FloatingPointError("Product receiver heat capacity became non-positive")
        temperature_rate = (
            enthalpy_rate_w - np.dot(component_energy, component_rate)
        ) / capacity
        return np.r_[component_rate, temperature_rate]


class Network:
    """One core for all schedules. Connections are simultaneous boundary fluxes."""

    def __init__(self, case, schedule):
        self.case = case
        self.schedule = schedule
        # Validate the physical topology before allocating state vectors.  In
        # particular, a 4-bed/8-step case must not silently become an 8-bed/
        # 8-step case (or vice versa) through a malformed serialized schedule.
        self.schedule_contract = validate_schedule_contract(schedule)
        self.species = case["species"]
        self.ns = len(self.species)
        self.nb = schedule["number_of_beds"]
        self.bed = Bed(case["bed"], self.species)
        self.bed_process_size = self.nb*self.bed.size
        self.receiver = None
        receiver_spec = case.get("product_receiver")
        if receiver_spec is not None:
            self.receiver = ProductReceiver(receiver_spec, self.bed)
        self.receiver_offset = self.bed_process_size
        self.receiver_size = 0 if self.receiver is None else self.receiver.size
        self.process_size = self.bed_process_size + self.receiver_size
        self.nledger = 3*self.ns + 5
        self.size = self.process_size + self.nledger
        self.feed = case["feed"]
        self.yfeed = np.asarray(self.feed["mole_fractions"], float)
        if self.yfeed.shape!=(self.ns,) or abs(self.yfeed.sum()-1) > 1e-9 or np.any(self.yfeed < 0):
            raise ValueError("Feed fractions must be nonnegative and sum to one")
        # The numerical engine is shared by independent gas systems.  H2 is
        # the historical default, but a source-specific binary CO2/N2
        # breakthrough case must be able to report CO2 without fabricating an
        # H2 component.  The declared reporting species remains the first
        # component so legacy product/recovery bookkeeping stays unchanged.
        self.reporting_species = case.get("reporting_species", "H2")
        names = [entry["name"] for entry in self.species]
        if self.reporting_species not in names:
            raise ValueError(f"Reporting species {self.reporting_species} is not present in the declared ordering")
        self.reporting_index = names.index(self.reporting_species)
        if self.case["low_pressure_bar"]<=0 or self.feed["pressure_bar"]<=self.case["low_pressure_bar"]:
            raise ValueError("Invalid high/low pressure ordering")
        self.valves = case["valves"]
        if case.get("boundary_mode") == "fixed_inlet_valve_outlet":
            product_valve = self.valves.get("product", {})
            if "area_cd_m2" not in product_valve:
                raise ValueError(
                    "fixed_inlet_valve_outlet requires an explicit product "
                    "orifice area_cd_m2; do not infer an outlet valve from flow"
                )
        if any(m.get("model") in {"streb_sips_wilson_rast_candidate", "guderian_2021_sips_iast"}
               for m in self.bed.materials):
            method = case.get("solver", {}).get("jacobian_method", "complex_step")
            if method not in {"bounded_fd", "scipy_fd"}:
                raise ValueError("Nonlinear source equilibrium requires bounded_fd or scipy_fd Jacobian")
        self.streb_topology_screen = schedule["name"] == "streb_two_bed_8operation_screen"
        self.streb_pressure_screen = case.get("boundary_mode") == "streb_pressure_equation_screen"
        self.guderian_topology = schedule["name"] == "guderian_two_bed_6step_source"
        if self.guderian_topology:
            if case.get("boundary_mode") != "guderian_receiver_cycle" or self.nb != 2:
                raise ValueError("Guderian source timing requires two beds and guderian_receiver_cycle boundaries")
            if self.receiver is None:
                raise ValueError("Guderian receiver cycle requires an explicit product_receiver state")
            boundary = case.get("source_boundary")
            if not isinstance(boundary, dict):
                raise ValueError("Guderian receiver cycle requires explicit source_boundary data")
            required_positive = (
                "purge_molar_flow_mol_s", "product_molar_flow_mol_s",
                "observed_cycle_mean_feed_mol_s", "ambient_tail_pressure_bar",
            )
            for key in required_positive:
                value = boundary.get(key)
                if value is None or not np.isfinite(value) or value <= 0.0:
                    raise ValueError(f"Invalid Guderian source boundary: {key}")
            if boundary.get("feed_boundary_mode") != "linear_pressure_driven_from_8.45_bar_supply":
                raise ValueError("Guderian source feed requires the pressure-driven Table-4 boundary")
            for name in ("feed", "tail", "receiver_connection", "equalization"):
                valve_spec = self.valves.get(name)
                if not isinstance(valve_spec, dict) or set(valve_spec) != {"gain_mol_s_bar"}:
                    raise ValueError(f"Guderian source valve {name} requires one linear gain_mol_s_bar")
                if not np.isfinite(valve_spec["gain_mol_s_bar"]) or valve_spec["gain_mol_s_bar"] <= 0.0:
                    raise ValueError(f"Guderian source valve {name} gain must be finite and positive")
            provenance = boundary.get("parameter_provenance")
            if not isinstance(provenance, str) or not provenance.strip():
                raise ValueError("Guderian source boundary requires parameter_provenance")
            self.source_boundary = dict(boundary)
            if self.bed.equation_of_state != "peng_robinson":
                raise ValueError("Guderian source adapter requires the published Peng-Robinson gas closure")
            if self.bed.energy_mode != "internal_energy":
                raise ValueError("Guderian receiver coupling requires conservative internal-energy mode")
        if self.streb_topology_screen:
            if case.get("boundary_mode") not in {"streb_topology_screen", "streb_pressure_equation_screen"} or self.nb != 2:
                raise ValueError("Streb timing requires two beds and explicit streb_topology_screen or streb_pressure_equation_screen boundaries")
            boundary = case.get("source_boundary")
            if not isinstance(boundary, dict):
                raise ValueError("Streb topology screen requires explicit source_boundary")
            required_positive = (
                "heavy_purge_molar_flow_mol_s", "light_purge_molar_flow_mol_s",
                "heavy_purge_temperature_k", "light_purge_temperature_k",
                "ambient_tail_pressure_bar",
            )
            for key in required_positive:
                value = boundary.get(key)
                if value is None or not np.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid Streb topology boundary: {key}")
            self.source_boundary = dict(boundary)
            heat_by_action = boundary.get("heat_transfer_w_m2_k_by_action")
            if heat_by_action is not None:
                expected_actions = {"PE-Press", "Feed", "Idle", "PE-BD", "BD", "HP", "BD-vac", "LP"}
                if (not isinstance(heat_by_action, dict)
                        or set(heat_by_action) != expected_actions):
                    raise ValueError("Streb step heat transfer requires all eight named actions")
                if any(not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0
                       for value in heat_by_action.values()):
                    raise ValueError("Streb step heat transfer coefficients must be positive and finite")
                provenance = boundary.get("heat_transfer_parameter_provenance")
                if not isinstance(provenance, str) or not provenance.strip():
                    raise ValueError("Streb step heat transfer requires heat_transfer_parameter_provenance")
                self.source_boundary["heat_transfer_w_m2_k_by_action"] = dict(heat_by_action)
            for label in ("heavy", "light"):
                key = f"{label}_purge_mole_fractions"
                composition = np.asarray(boundary.get(key), dtype=float)
                if (composition.shape != (self.ns,) or not np.all(np.isfinite(composition))
                        or np.any(composition < 0) or abs(composition.sum() - 1) > 1e-9):
                    raise ValueError(f"Invalid Streb topology boundary: {key}")
                self.source_boundary[key] = composition
            if self.streb_pressure_screen:
                if self.bed.mode != "loading_ldf":
                    raise ValueError("Streb pressure-equation screen requires loading_ldf kinetics")
                if case.get("solver", {}).get("jacobian_method", "complex_step") not in {"bounded_fd", "scipy_fd"}:
                    raise ValueError("Streb pressure-equation screen requires bounded_fd or scipy_fd Jacobian")
                for action in ("bd", "bd_vac"):
                    for suffix in ("asymptote_bar", "decay_per_s"):
                        key = f"{action}_{suffix}"
                        value = boundary.get(key)
                        if value is None or not np.isfinite(value) or value <= 0:
                            raise ValueError(f"Invalid Streb pressure-equation boundary: {key}")
                if self.source_boundary["bd_vac_asymptote_bar"] >= self.source_boundary["bd_asymptote_bar"]:
                    raise ValueError("BD-vac asymptote must be below BD asymptote")
                if abs(self.source_boundary["bd_asymptote_bar"] - self.source_boundary["ambient_tail_pressure_bar"]) > 1e-9:
                    raise ValueError("BD asymptote must match the declared ambient tail pressure")
                if not isinstance(boundary.get("pressure_parameter_provenance"), str) or not boundary["pressure_parameter_provenance"].strip():
                    raise ValueError("Streb pressure-equation boundary requires pressure_parameter_provenance")
                lp_keys = {"lp_initial_pressure_bar", "lp_final_pressure_bar",
                           "lp_pressure_parameter_provenance"}
                self.streb_lp_pressure_screen = bool(lp_keys & boundary.keys())
                if self.streb_lp_pressure_screen:
                    if not lp_keys <= boundary.keys():
                        raise ValueError("Streb LP linear pressure requires initial, final and provenance")
                    start = boundary["lp_initial_pressure_bar"]
                    end = boundary["lp_final_pressure_bar"]
                    if (not isinstance(start, (int, float)) or not isinstance(end, (int, float))
                            or not np.isfinite(start) or not np.isfinite(end)
                            or not 0 < end < start):
                        raise ValueError("Streb LP linear pressure requires initial > final > 0")
                    if (not isinstance(boundary["lp_pressure_parameter_provenance"], str)
                            or not boundary["lp_pressure_parameter_provenance"].strip()):
                        raise ValueError("Streb LP linear pressure requires nonempty provenance")
                    lp_durations = [sum(slot["duration_s"] for slot in schedule["slots"]
                                        if slot["actions"][bed] == "LP") for bed in range(self.nb)]
                    if (not all(value > 0 for value in lp_durations)
                            or abs(lp_durations[0]-lp_durations[1]) > 1e-9):
                        raise ValueError("Streb LP duration must be positive and equal in both beds")
                    self.streb_lp_duration_s = lp_durations[0]
                else:
                    self.streb_lp_duration_s = None
        elif (case.get("boundary_mode") in {"streb_topology_screen", "streb_pressure_equation_screen"}
              or "heat_transfer_w_m2_k_by_action" in case.get("source_boundary", {})):
            raise ValueError("Streb topology boundaries require the dedicated two-bed schedule")
        if (case.get("boundary_mode") == "guderian_receiver_cycle") != self.guderian_topology:
            raise ValueError("Guderian receiver boundaries require the dedicated source schedule")

    def initial(self):
        y = np.asarray(self.case.get("initial_mole_fractions", []), dtype=float)
        if y.size == 0:
            y = np.zeros(self.ns); y[0] = 1.0
        if y.shape != (self.ns,) or not np.all(np.isfinite(y)) or np.any(y < 0) or not np.isclose(y.sum(), 1.0, atol=1e-10):
            raise ValueError("Initial mole fractions must be nonnegative, finite, and sum to one")
        pressures = self.case.get("initial_bed_pressures_bar")
        compositions = self.case.get("initial_bed_mole_fractions")
        if pressures is None:
            pressures = [self.case["initial_pressure_bar"]] * self.nb
        if compositions is None:
            compositions = [y.tolist()] * self.nb
        if len(pressures) != self.nb or len(compositions) != self.nb:
            raise ValueError("Initial bed pressures and compositions must have one entry per physical bed")
        beds = []
        for index, (pressure, composition) in enumerate(zip(pressures, compositions)):
            composition = np.asarray(composition, dtype=float)
            if (not np.isfinite(pressure) or pressure <= 0.0
                    or composition.shape != (self.ns,) or np.any(~np.isfinite(composition))
                    or np.any(composition < 0.0) or not np.isclose(composition.sum(), 1.0, atol=1e-10)):
                raise ValueError(f"Invalid initial state for bed {index}")
            beds.append(self.bed.initial(
                float(pressure), composition, self.case["bed"]["ambient_k"]
            ))
        receiver = []
        if self.receiver is not None:
            receiver = self.receiver.initial(
                self.case["product_receiver"]["initial_pressure_bar"],
                self.case["product_receiver"]["initial_mole_fractions"],
                self.case["product_receiver"]["initial_temperature_k"],
            )
        return np.r_[np.concatenate(beds), receiver, np.zeros(self.nledger)]

    def inventory_energy(self, state):
        inv = np.zeros(self.ns); energy = 0.
        for i in range(self.nb):
            n, e = self.bed.inventory_energy(state[i*self.bed.size:(i+1)*self.bed.size])
            inv += n; energy += e
        if self.receiver is not None:
            n, e = self.receiver.inventory_energy(
                state[self.receiver_offset:self.receiver_offset+self.receiver_size]
            )
            inv += n; energy += e
        return inv, energy

    def rhs(self, time, state, slot, ledger_only=False, diagnostics=None):
        """Evaluate the process RHS and boundary ledger.

        ``diagnostics`` is an optional read-only audit sink.  It is deliberately
        opt-in so normal solver calls, including complex-step Jacobian calls,
        retain exactly the historical behavior.  When supplied, scalar
        boundary-flow contributions are recorded by bed/action for diagnosing
        source-apparatus boundary semantics; they do not feed back into the
        state equations.
        """
        n = self.bed.n
        states = [state[i*self.bed.size:(i+1)*self.bed.size] for i in range(self.nb)]
        fields = [self.bed.fields(x) for x in states]
        receiver_state = None
        receiver_fields = None
        receiver_component_rate = None
        receiver_enthalpy_rate = 0.0
        if self.receiver is not None:
            receiver_state = state[self.receiver_offset:self.receiver_offset+self.receiver_size]
            receiver_fields = self.receiver.fields(receiver_state)
            receiver_component_rate = np.zeros(self.ns, dtype=np.result_type(state, float))
        dtype=np.result_type(state,float)
        fluxes = ([(np.zeros((n+1,self.ns),dtype=dtype),np.zeros(n+1,dtype=dtype)) for _ in states]
                  if ledger_only and not self.streb_pressure_screen
                  else [self.bed.internal_flux(x) for x in states])
        ledger = np.zeros(self.nledger,dtype=dtype)
        if diagnostics is not None:
            diagnostics.clear()
            diagnostics.update({
                "feed_flow_by_bed_mol_s": [0.0]*self.nb,
                "adsorption_gross_flow_by_bed_mol_s": [0.0]*self.nb,
                "tail_flow_by_bed_mol_s": [0.0]*self.nb,
                "tail_flow_by_action_mol_s": {},
                "purge_received_flow_by_bed_mol_s": [0.0]*self.nb,
                "purge_donated_flow_by_bed_mol_s": [0.0]*self.nb,
                "product_exported_flow_by_bed_mol_s": [0.0]*self.nb,
                "equalization_flow_by_pair_mol_s": [],
                "product_purge_flow_by_pair_mol_s": [],
                "purge_requested_mol_s": 0.0,
                "purge_unmet_mol_s": 0.0,
                "source_external_inflow_by_action_mol_s": {},
                "source_external_outflow_by_action_mol_s": {},
                "source_pressure_control_by_bed": [],
                "receiver_pressure_bar": (None if receiver_fields is None
                                          else float(np.real(receiver_fields[2]))),
                "receiver_temperature_k": (None if receiver_fields is None
                                            else float(np.real(receiver_fields[1]))),
                "receiver_flow_events": [],
                "receiver_product_export_mol_s": 0.0,
                "receiver_purge_requested_mol_s": 0.0,
                "receiver_purge_delivered_mol_s": 0.0,
                "valve_resolution": {},
            })
        actions = slot["actions"]
        heat_by_action = (self.source_boundary.get("heat_transfer_w_m2_k_by_action")
                          if self.streb_topology_screen else None)

        def gas_wall_coefficient(bed_index):
            return (self.bed.spec.get("h_bed_wall_w_m2_k", 0.0)
                    if heat_by_action is None else heat_by_action[actions[bed_index]])
        ad_streams = []
        controlled_outflows = []

        def apply(bed, endpoint, flow, yy, temperature, pressure, incoming):
            face = 0 if endpoint == "bottom" else n
            signed = flow * (1 if incoming else -1) * (1 if face == 0 else -1)
            fluxes[bed][0][face] += signed*yy
            fluxes[bed][1][face] += signed*self.bed.enthalpy(yy, temperature, pressure)

        def valve(pu, pd, yy, tt, name):
            # Product-purge is a distinct source operation in the Linde
            # eight-step topology.  A case must provide its own valve law;
            # falling back to the equalization law is explicit and only keeps
            # legacy exploratory cases executable.
            resolved_name = name if name in self.valves else "equalization"
            spec = self.valves.get(resolved_name)
            if spec is None:
                raise KeyError(f"Missing valve specification: {name}")
            if diagnostics is not None:
                diagnostics["valve_resolution"][name] = {
                    "resolved_as": resolved_name,
                    "fallback_used": resolved_name != name,
                }
            regularization=self.case.get("valve_regularization_bar",1e-6)
            if "area_cd_m2" in spec:
                return valve_molar_flow(pu, pd, tt, yy, self.bed.mw, self.bed.cp, spec["area_cd_m2"],regularization)
            delta=_positive(pu-pd)
            if "gain_mol_s_bar" in spec:
                return spec["gain_mol_s_bar"]*delta
            return spec["gain_mol_s_sqrt_bar"]*delta/np.sqrt(delta if np.real(delta)>regularization else regularization)

        def source_in(bed, action, endpoint, flow, yy, temperature, pressure=None):
            if pressure is None:
                pressure = self.feed["pressure_bar"]
            apply(bed, endpoint, flow, yy, temperature, pressure, True)
            # The legacy feed ledger represents every external inlet.  The
            # diagnostic action tag keeps the feed/HP/LP origins distinct.
            ledger[:self.ns] += flow*yy
            ledger[3*self.ns] += flow*self.bed.enthalpy(yy, temperature, pressure)
            if diagnostics is not None:
                flows = diagnostics["source_external_inflow_by_action_mol_s"]
                flows[action] = flows.get(action, 0.0) + float(np.real(flow))
                if action in {"Feed", "AP"}:
                    diagnostics["feed_flow_by_bed_mol_s"][bed] += float(np.real(flow))

        def source_out(bed, action, endpoint, pressure, yy, temperature,
                       product=False, forced_flow=None):
            if product:
                target = self.case["product_pressure_bar"]
            elif action in {"BD", "HP"}:
                target = self.source_boundary["ambient_tail_pressure_bar"]
            elif self.streb_pressure_screen and action == "BD-vac":
                target = self.source_boundary["bd_vac_asymptote_bar"]
            elif self.streb_pressure_screen and action == "LP" and self.streb_lp_pressure_screen:
                target = self.source_boundary["lp_final_pressure_bar"]
            else:
                target = self.case["low_pressure_bar"]
            flow = (forced_flow if forced_flow is not None else
                    valve(pressure, target, yy, temperature, "product" if product else "tail"))
            apply(bed, endpoint, flow, yy, temperature, pressure, False)
            if product:
                ledger[self.ns:2*self.ns] += flow*yy
                ledger[3*self.ns+1] += flow*self.bed.enthalpy(yy, temperature, pressure)
            else:
                ledger[2*self.ns:3*self.ns] += flow*yy
                ledger[3*self.ns+2] += flow*self.bed.enthalpy(yy, temperature, pressure)
                if target < 1.01325:
                    ledger[-1] += (flow*R*temperature*math.log(1.01325/target)
                                   /self.case.get("vacuum_efficiency", 0.65))
            if diagnostics is not None:
                flows = diagnostics["source_external_outflow_by_action_mol_s"]
                flows[action] = flows.get(action, 0.0) + float(np.real(flow))
                key = ("product_exported_flow_by_bed_mol_s" if product
                       else "tail_flow_by_bed_mol_s")
                diagnostics[key][bed] += float(np.real(flow))
            return flow

        def receiver_transfer_to_bed(bed, endpoint, flow, action):
            nonlocal receiver_enthalpy_rate
            if self.receiver is None:
                raise RuntimeError("Receiver transfer requested without a receiver state")
            _, temperature, pressure, composition = receiver_fields
            apply(bed, endpoint, flow, composition, temperature, pressure, True)
            receiver_component_rate[:] -= flow * composition
            receiver_enthalpy_rate -= flow * self.receiver.enthalpy(composition, temperature, pressure)
            if diagnostics is not None:
                diagnostics["receiver_flow_events"].append({
                    "action": action, "direction": "receiver_to_bed", "bed": int(bed),
                    "endpoint": endpoint, "flow_mol_s": float(np.real(flow)),
                })

        def bed_transfer_to_receiver(bed, endpoint, flow, composition, temperature, pressure, action):
            nonlocal receiver_enthalpy_rate
            if self.receiver is None:
                raise RuntimeError("Receiver transfer requested without a receiver state")
            apply(bed, endpoint, flow, composition, temperature, pressure, False)
            receiver_component_rate[:] += flow * composition
            receiver_enthalpy_rate += flow * self.bed.enthalpy(composition, temperature, pressure)
            if diagnostics is not None:
                diagnostics["receiver_flow_events"].append({
                    "action": action, "direction": "bed_to_receiver", "bed": int(bed),
                    "endpoint": endpoint, "flow_mol_s": float(np.real(flow)),
                })

        def bidirectional_receiver_connection(bed, endpoint, pressure, composition,
                                              temperature, action):
            receiver_pressure = receiver_fields[2]
            if np.real(receiver_pressure) >= np.real(pressure):
                flow = valve(receiver_pressure, pressure, receiver_fields[3],
                             receiver_fields[1], "receiver_connection")
                receiver_transfer_to_bed(bed, endpoint, flow, action)
                direction = "receiver_to_bed"
            else:
                flow = valve(pressure, receiver_pressure, composition, temperature,
                             "receiver_connection")
                bed_transfer_to_receiver(bed, endpoint, flow, composition, temperature, pressure, action)
                direction = "bed_to_receiver"
            return flow, direction

        def receiver_product_export():
            nonlocal receiver_enthalpy_rate
            if self.receiver is None:
                return
            inventory, temperature, pressure, composition = receiver_fields
            requested = float(self.source_boundary["product_molar_flow_mol_s"])
            # A source flow controller cannot withdraw more material than is
            # present.  The small inventory floor is a numerical safeguard;
            # an unmet demand remains visible in the trajectory and prevents
            # source-equivalent qualification.
            available = float(np.real(inventory.sum())) / max(
                float(self.schedule["cycle_time_s"]), 1.0
            )
            flow = min(requested, max(available, 0.0))
            receiver_component_rate[:] -= flow * composition
            enthalpy = self.receiver.enthalpy(composition, temperature, pressure)
            receiver_enthalpy_rate -= flow * enthalpy
            ledger[self.ns:2*self.ns] += flow * composition
            ledger[3*self.ns+1] += flow * enthalpy
            if diagnostics is not None:
                diagnostics["receiver_product_export_mol_s"] = float(flow)
                diagnostics["product_exported_flow_by_bed_mol_s"][0] += 0.0

        for i, action in enumerate(actions):
            _, _, t, _, p, y = fields[i]
            if self.guderian_topology:
                if action == "AP":
                    # The source feed is a compressed-air pressure boundary
                    # followed by the V1/V2(or V3)/support/microfilter linear
                    # resistance chain.  The workbook flow is a measured
                    # cycle-average target, not an instantaneous imposed MFC.
                    feed_flow = valve(
                        self.feed["pressure_bar"], p[0], self.yfeed,
                        self.feed["temperature_k"], "feed"
                    )
                    source_in(i, action, "bottom", feed_flow,
                              self.yfeed, self.feed["temperature_k"])
                    bidirectional_receiver_connection(
                        i, "top", p[-1], y[-1], t[-1], action
                    )
                elif action == "BR":
                    source_out(i, action, "bottom", p[0], y[0], t[0])
                    requested = float(self.source_boundary["purge_molar_flow_mol_s"])
                    receiver_pressure = receiver_fields[2]
                    delivered = requested if np.real(receiver_pressure) > np.real(p[-1]) else 0.0
                    receiver_transfer_to_bed(i, "top", delivered, action)
                    if diagnostics is not None:
                        diagnostics["receiver_purge_requested_mol_s"] += requested
                        diagnostics["receiver_purge_delivered_mol_s"] += delivered
                        diagnostics["purge_requested_mol_s"] += requested
                        diagnostics["purge_unmet_mol_s"] += requested - delivered
                        diagnostics["purge_received_flow_by_bed_mol_s"][i] += delivered
                elif action not in {"EQ-D", "EQ-R"}:
                    raise ValueError(f"Unsupported Guderian source action: {action}")
                continue
            if self.streb_topology_screen:
                # Source topology only.  Gas-bottle MFC inlets and valve
                # outflows do not reproduce the paper's measured pressure
                # trajectories or its Wilson/RAST equilibrium closure.
                if action == "Feed":
                    source_in(i, action, "bottom", self.feed["molar_flow_mol_s"],
                              self.yfeed, self.feed["temperature_k"])
                    source_out(i, action, "top", p[-1], y[-1], t[-1], product=True)
                elif action == "BD":
                    if self.streb_pressure_screen:
                        controlled_outflows.append((i, action, "top"))
                    else:
                        source_out(i, action, "top", p[-1], y[-1], t[-1])
                elif action == "HP":
                    source_in(i, action, "bottom",
                              self.source_boundary["heavy_purge_molar_flow_mol_s"],
                              self.source_boundary["heavy_purge_mole_fractions"],
                              self.source_boundary["heavy_purge_temperature_k"])
                    source_out(i, action, "top", p[-1], y[-1], t[-1])
                elif action == "BD-vac":
                    if self.streb_pressure_screen:
                        controlled_outflows.append((i, action, "bottom"))
                    else:
                        source_out(i, action, "bottom", p[0], y[0], t[0])
                elif action == "LP":
                    source_in(i, action, "top",
                              self.source_boundary["light_purge_molar_flow_mol_s"],
                              self.source_boundary["light_purge_mole_fractions"],
                              self.source_boundary["light_purge_temperature_k"])
                    if self.streb_pressure_screen and self.streb_lp_pressure_screen:
                        controlled_outflows.append((i, action, "bottom"))
                    else:
                        source_out(i, action, "bottom", p[0], y[0], t[0])
                elif action not in {"PE-Press", "PE-BD", "Idle"}:
                    raise ValueError(f"Unsupported Streb topology action: {action}")
                continue
            if action in {"FP", "AD"}:
                spec = self.valves["feed"]
                if self.case.get("boundary_mode") in {
                    "fixed_flow_breakthrough", "fixed_inlet_valve_outlet"
                }:
                    # A fixed-flow breakthrough test uses the reported MFC
                    # rate as an explicit boundary condition.  This is not a
                    # plant valve rating; the case must separately document
                    # its pressure-control and outlet assumptions.
                    flow = self.feed["molar_flow_mol_s"]
                elif self.case.get("boundary_mode") == "isobaric_breakthrough":
                    # A laboratory breakthrough bed is held near its
                    # reported pressure by a back-pressure regulator.  The
                    # finite-volume state still carries adsorption uptake,
                    # so a bounded proportional inlet correction supplies
                    # the missing gas without replacing the transport ODE.
                    gain = self.case.get("pressure_control_gain_mol_s_bar", 0.01)
                    flow = self.feed["molar_flow_mol_s"] + gain * _positive(self.feed["pressure_bar"] - p[0])
                elif self.case.get("boundary_mode") == "wang_controller":
                    taper = _minimum(1.,_positive((self.feed["pressure_bar"]-p[0])/spec["taper_band_bar"])) if action == "FP" else 1.
                    flow = self.feed["molar_flow_mol_s"] * taper
                else:
                    flow = _minimum(self.feed["molar_flow_mol_s"], valve(self.feed["pressure_bar"], p[0], self.yfeed, self.feed["temperature_k"], "feed"))
                apply(i, "bottom", flow, self.yfeed, self.feed["temperature_k"], self.feed["pressure_bar"], True)
                ledger[:self.ns] += flow*self.yfeed
                ledger[3*self.ns] += flow*self.bed.enthalpy(self.yfeed, self.feed["temperature_k"], self.feed["pressure_bar"])
                if diagnostics is not None:
                    diagnostics["feed_flow_by_bed_mol_s"][i] += float(np.real(flow))
            if action == "AD":
                if self.case.get("boundary_mode") == "fixed_inlet_valve_outlet":
                    # Fixed MFC inlet with a declared physical outlet valve.
                    # This is a generic boundary option for source protocols
                    # that report an inlet MFC but not a prescribed outlet
                    # flow.  It is never source-equivalent without the source
                    # valve/boundary data.
                    out = valve(p[-1], self.case["product_pressure_bar"],
                                y[-1], t[-1], "product")
                elif self.case.get("boundary_mode") in {"fixed_flow_breakthrough", "isobaric_breakthrough"}:
                    out = self.feed["molar_flow_mol_s"]
                elif self.case.get("boundary_mode") == "wang_controller":
                    if self.case.get("wang_product_boundary") == "fixed_feed_flow_diagnostic":
                        # Read-only source-boundary diagnostic: the Wang
                        # apparatus reports a fixed MFC feed rate.  Keep the
                        # historical pressure-controller law as the default;
                        # this branch is never used by frozen cases.
                        out = self.feed["molar_flow_mol_s"]
                    else:
                        out = _positive(self.feed["molar_flow_mol_s"] + self.valves["product"]["gain_mol_s_bar"]*(p[-1]-self.case["high_pressure_bar"]))
                else:
                    out = valve(p[-1], self.case["product_pressure_bar"], y[-1], t[-1], "product")
                ad_streams.append([i, out, y[-1], t[-1], p[-1]])
                if diagnostics is not None:
                    diagnostics["adsorption_gross_flow_by_bed_mol_s"][i] += float(np.real(out))
            if action in {"BD", "P", "EVC"}:
                if self.case.get("boundary_mode") == "wang_controller":
                    out = self.valves["tail"]["gain_mol_s_bar"]*_positive(p[0]-self.case["low_pressure_bar"])
                else:
                    out = valve(p[0], self.case["low_pressure_bar"], y[0], t[0], "tail")
                apply(i, "bottom", out, y[0], t[0], p[0], False)
                ledger[2*self.ns:3*self.ns] += out*y[0]
                ledger[3*self.ns+2] += out*self.bed.enthalpy(y[0], t[0], p[0])
                if diagnostics is not None:
                    diagnostics["tail_flow_by_bed_mol_s"][i] += float(np.real(out))
                    by_action = diagnostics["tail_flow_by_action_mol_s"]
                    by_action[action] = by_action.get(action, 0.0) + float(np.real(out))
                if self.case["low_pressure_bar"] < 1.01325:
                    ledger[-1] += out*R*t[0]*math.log(1.01325/self.case["low_pressure_bar"])/self.case.get("vacuum_efficiency", 0.65)

        # Purge is taken from the instantaneous real adsorption outlet. It is
        # subtracted before external product is counted, including startup.
        purge_beds = [i for i, a in enumerate(actions) if a == "P"]
        for target in purge_beds:
            remaining = self.feed["molar_flow_mol_s"]*self.case["purge_to_feed_ratio"]
            requested = float(np.real(remaining))
            for stream in ad_streams:
                donor, available, yy, tt, pp = stream
                target_p = fields[target][4][-1]
                donor_p = fields[donor][4][-1]
                capacity = _positive(donor_p-target_p)*self.case.get("purge_gain_mol_s_bar", 0.01)
                flow = _minimum(remaining, available, capacity)
                apply(target, "top", flow, yy, tt, pp, True)
                stream[1] -= flow
                remaining -= flow
                if diagnostics is not None:
                    diagnostics["purge_received_flow_by_bed_mol_s"][target] += float(np.real(flow))
                    diagnostics["purge_donated_flow_by_bed_mol_s"][donor] += float(np.real(flow))
            if diagnostics is not None:
                diagnostics["purge_requested_mol_s"] += requested
                diagnostics["purge_unmet_mol_s"] += float(np.real(_positive(remaining)))

        for donor, exported, yy, tt, pp in ad_streams:
            # Restore the gross outflow: exported plus all transfers from this
            # donor. Its top flux must equal the product + purge inventory.
            if self.case.get("boundary_mode") == "fixed_inlet_valve_outlet":
                gross = valve(fields[donor][4][-1], self.case["product_pressure_bar"],
                              yy, tt, "product")
            elif self.case.get("boundary_mode") in {"fixed_flow_breakthrough", "isobaric_breakthrough"}:
                gross = self.feed["molar_flow_mol_s"]
            elif self.case.get("boundary_mode") == "wang_controller":
                if self.case.get("wang_product_boundary") == "fixed_feed_flow_diagnostic":
                    gross = self.feed["molar_flow_mol_s"]
                else:
                    gross = _positive(self.feed["molar_flow_mol_s"] + self.valves["product"]["gain_mol_s_bar"]*(fields[donor][4][-1]-self.case["high_pressure_bar"]))
            else:
                gross = valve(fields[donor][4][-1], self.case["product_pressure_bar"], yy, tt, "product")
            apply(donor, "top", gross, yy, tt, pp, False)
            ledger[self.ns:2*self.ns] += exported*yy
            ledger[3*self.ns+1] += exported*self.bed.enthalpy(yy, tt, pp)
            if diagnostics is not None:
                diagnostics["product_exported_flow_by_bed_mol_s"][donor] += float(np.real(exported))

        # The source Linde table distinguishes inter-bed product purge from
        # the historical P action.  PPD/PPR are explicit boundary states, so
        # this transfer is not counted as external product and is never
        # silently replaced by the automatic adsorption-outlet purge logic.
        for donor, receiver in slot.get("product_purge_pairs", []):
            pd, pr = fields[donor][4][-1], fields[receiver][4][-1]
            yy, tt = fields[donor][5][-1], fields[donor][2][-1]
            flow = valve(pd, pr, yy, tt, "product_purge")
            apply(donor, "top", flow, yy, tt, pd, False)
            apply(receiver, "top", flow, yy, tt, pd, True)
            if diagnostics is not None:
                diagnostics["product_purge_flow_by_pair_mol_s"].append({
                    "donor": int(donor), "receiver": int(receiver),
                    "flow_mol_s": float(np.real(flow))})

        for donor, receiver in slot.get("equalization_pairs", []):
            pd, pr = fields[donor][4][-1], fields[receiver][4][-1]
            yy, tt = fields[donor][5][-1], fields[donor][2][-1]
            flow = valve(pd, pr, yy, tt, "equalization")
            apply(donor, "top", flow, yy, tt, pd, False)
            apply(receiver, "top", flow, yy, tt, pd, True)
            if diagnostics is not None:
                diagnostics["equalization_flow_by_pair_mol_s"].append({
                    "donor": int(donor), "receiver": int(receiver),
                    "flow_mol_s": float(np.real(flow)),
                    # The same donor product-end composition used by the
                    # physical inter-bed flux. Diagnostic-only: this does
                    # not add an external stream or change either bed RHS.
                    "species_flow_mol_s": np.real(flow * yy).astype(float).tolist()})

        # The Guderian apparatus opens both like-end equalisation paths at the
        # same physical time.  Direction follows the instantaneous pressure,
        # even if a numerical perturbation reverses the nominal donor label.
        for connection in slot.get("equalization_connections", []):
            first = int(connection["donor"])
            second = int(connection["receiver"])
            first_endpoint = connection["donor_endpoint"]
            second_endpoint = connection["receiver_endpoint"]
            first_cell = 0 if first_endpoint == "bottom" else n - 1
            second_cell = 0 if second_endpoint == "bottom" else n - 1
            first_pressure = fields[first][4][first_cell]
            second_pressure = fields[second][4][second_cell]
            if np.real(first_pressure) >= np.real(second_pressure):
                actual_donor, actual_receiver = first, second
                donor_endpoint, receiver_endpoint = first_endpoint, second_endpoint
                donor_cell = first_cell
            else:
                actual_donor, actual_receiver = second, first
                donor_endpoint, receiver_endpoint = second_endpoint, first_endpoint
                donor_cell = second_cell
            pressure_up = fields[actual_donor][4][donor_cell]
            receiver_cell = 0 if receiver_endpoint == "bottom" else n - 1
            pressure_down = fields[actual_receiver][4][receiver_cell]
            composition = fields[actual_donor][5][donor_cell]
            temperature = fields[actual_donor][2][donor_cell]
            flow = valve(pressure_up, pressure_down, composition, temperature, "equalization")
            apply(actual_donor, donor_endpoint, flow, composition, temperature, pressure_up, False)
            apply(actual_receiver, receiver_endpoint, flow, composition, temperature, pressure_up, True)
            if diagnostics is not None:
                diagnostics["equalization_flow_by_pair_mol_s"].append({
                    "donor": int(actual_donor), "receiver": int(actual_receiver),
                    "donor_endpoint": donor_endpoint,
                    "receiver_endpoint": receiver_endpoint,
                    "flow_mol_s": float(np.real(flow)),
                    "species_flow_mol_s": np.real(flow * composition).astype(float).tolist(),
                })

        if self.guderian_topology:
            receiver_product_export()

        # In the source-inspired diagnostic, prescribe the *derivative* of
        # endpoint-cell-center pressure during BD/BD-vac and optional LP, not
        # the gas inventory
        # itself. This is not the paper's Dirichlet pressure at the port face.
        # The required external outflow follows from the existing finite-
        # volume mass/energy RHS.  It is booked through source_out exactly
        # like a physical boundary flow, preserving both conservation ledgers.
        # No sign reversal, hidden makeup gas or guessed valve rating is used.
        for bed_index, action, endpoint in controlled_outflows:
            cell = 0 if endpoint == "bottom" else n - 1
            gas, _, temperatures, _, pressures, fractions = fields[bed_index]
            pressure = float(np.real(pressures[cell]))
            temperature = float(np.real(temperatures[cell]))
            if action == "LP":
                from .streb_2021_pressure import linear_pressure_bar
                target_slope = linear_pressure_bar(
                    self.source_boundary["lp_initial_pressure_bar"],
                    self.source_boundary["lp_final_pressure_bar"],
                    self.streb_lp_duration_s, 0.0)[1]
            else:
                suffix = "bd" if action == "BD" else "bd_vac"
                asymptote = float(self.source_boundary[f"{suffix}_asymptote_bar"])
                decay = float(self.source_boundary[f"{suffix}_decay_per_s"])
                if pressure <= asymptote:
                    raise ValueError(f"{action} pressure must exceed its declared asymptote")
                target_slope = exponential_pressure_bar(pressure, asymptote, decay, 0.0)[1]
            baseline, _ = self.bed.rates(states[bed_index], *fluxes[bed_index],
                                         h_bed_wall_override=gas_wall_coefficient(bed_index))
            baseline = baseline.reshape(n, self.bed.nv)
            ngas = float(np.real(np.sum(gas[cell])))
            base_slope = R / (self.bed.vgas[cell] * BAR) * (
                temperature * float(np.real(np.sum(baseline[cell, :self.ns])))
                + ngas * float(np.real(baseline[cell, -2])))
            composition = fractions[cell]
            outgoing_enthalpy = float(np.real(self.bed.enthalpy(composition, temperature, pressure)))
            capacity = float(self.bed.cs[cell])
            numerator = -outgoing_enthalpy
            if self.bed.energy_mode == "internal_energy":
                capacity += float(np.real(np.dot(gas[cell], self.bed.cv)))
                numerator += temperature * float(np.real(np.dot(composition, self.bed.cv)))
            unit_temperature_slope = numerator / capacity
            unit_pressure_slope = R / (self.bed.vgas[cell] * BAR) * (
                -temperature + ngas * unit_temperature_slope)
            if not np.isfinite(unit_pressure_slope) or unit_pressure_slope >= 0:
                raise ArithmeticError(f"{action} outflow does not lower pressure")
            flow = (target_slope - base_slope) / unit_pressure_slope
            if not np.isfinite(flow) or flow < -1e-12:
                raise ValueError(f"{action} pressure trajectory would require an inlet, not an outflow")
            flow = max(0.0, flow)
            source_out(bed_index, action, endpoint, pressure, composition,
                       temperature, forced_flow=flow)
            if diagnostics is not None:
                diagnostics["source_pressure_control_by_bed"].append({
                    "bed": bed_index, "action": action, "endpoint": endpoint,
                    "controlled_location": "endpoint_cell_center_not_port_face",
                    "pressure_bar": pressure, "target_dpdt_bar_s": target_slope,
                    "baseline_dpdt_bar_s": base_slope,
                    "outflow_mol_s": flow,
                    "predicted_dpdt_bar_s": base_slope + flow * unit_pressure_slope,
                })

        if ledger_only:
            if self.bed.wall_model == "fixed_ambient":
                ledger[3*self.ns+3] = sum(
                    gas_wall_coefficient(index) * np.sum(
                        self.bed.wall_area * (self.bed.spec["ambient_k"]-field[2])
                    )
                    for index, field in enumerate(fields)
                )
            else:
                ledger[3*self.ns+3]=sum(self.bed.spec.get("h_ambient_w_m2_k",0)*
                    np.sum(self.bed.outer_area*(self.bed.spec["ambient_k"]-field[3])) for field in fields)
            return ledger
        rates = []
        for i in range(self.nb):
            dx, ambient = self.bed.rates(states[i], *fluxes[i],
                                         h_bed_wall_override=gas_wall_coefficient(i))
            rates.append(dx)
            ledger[3*self.ns+3] += ambient
        receiver_rate = []
        if self.receiver is not None:
            receiver_rate = self.receiver.rates(
                receiver_state, receiver_component_rate, receiver_enthalpy_rate
            )
        return np.r_[np.concatenate(rates), receiver_rate, ledger]

    def sparsity(self,slot=None):
        pattern = lil_matrix((self.size, self.size), dtype=bool)
        nv = self.bed.nv
        for bed in range(self.nb):
            offset = bed*self.bed.size
            for cell in range(self.bed.n):
                rows = list(range(offset+cell*nv, offset+(cell+1)*nv))
                cols = list(range(offset+max(0, cell-1)*nv, offset+min(self.bed.n, cell+2)*nv))
                pattern[np.ix_(rows, cols)] = True
        boundaries = []
        for bed in range(self.nb):
            base = bed*self.bed.size
            boundaries += list(range(base, base+nv)) + list(range(base+self.bed.size-nv, base+self.bed.size))
        receiver_columns = list(range(self.receiver_offset, self.receiver_offset+self.receiver_size))
        if slot is None:
            pattern[np.ix_(boundaries, boundaries)] = True
        else:
            # Only physically connected top ports couple separate beds.
            # A global boundary clique unnecessarily gives every endpoint a
            # different derivative color even while those valves are closed.
            top=lambda bed:list(range((bed+1)*self.bed.size-nv,(bed+1)*self.bed.size))
            bottom=lambda bed:list(range(bed*self.bed.size,bed*self.bed.size+nv))
            for donor,receiver in slot.get("equalization_pairs",[]):
                connected=top(donor)+top(receiver)
                pattern[np.ix_(connected,connected)]=True
            for donor, receiver in slot.get("product_purge_pairs", []):
                connected=top(donor)+top(receiver)
                pattern[np.ix_(connected,connected)]=True
            for connection in slot.get("equalization_connections", []):
                donor = int(connection["donor"])
                receiver = int(connection["receiver"])
                donor_columns = (bottom(donor) if connection["donor_endpoint"] == "bottom"
                                 else top(donor))
                receiver_endpoint_columns = (
                    bottom(receiver) if connection["receiver_endpoint"] == "bottom"
                    else top(receiver)
                )
                connected = donor_columns + receiver_endpoint_columns
                pattern[np.ix_(connected, connected)] = True
            if self.receiver is not None:
                connected = list(receiver_columns)
                for bed, action in enumerate(slot["actions"]):
                    if action in {"AP", "BR"}:
                        connected += top(bed)
                pattern[np.ix_(connected, connected)] = True
            purge=[i for i,a in enumerate(slot["actions"]) if a=="P"]
            adsorption=[i for i,a in enumerate(slot["actions"]) if a=="AD"]
            if purge:
                connected=[col for bed in purge+adsorption for col in top(bed)]
                pattern[np.ix_(connected,connected)]=True
        pattern[self.process_size:, boundaries + receiver_columns] = True
        # Only wall temperatures enter the ambient heat ledger. A false dense
        # lower block needlessly makes every bed-state a different FD color.
        wall_columns = [b*self.bed.size+c*nv+nv-1 for b in range(self.nb) for c in range(self.bed.n)]
        pattern[self.process_size+3*self.ns+3, wall_columns] = True
        return pattern.tocsr()
