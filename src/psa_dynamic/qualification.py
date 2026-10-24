"""Fail-closed engineering evidence screens, distinct from ODE convergence.

A source DOI alone does not validate an extrapolated material or a scaled-up
vessel. These reports expose missing evidence without prohibiting clearly
labelled development simulations. They never authorize data release.
"""
from __future__ import annotations
import math
import numpy as np


class MaterialDomainError(ValueError):
    """Constitutive-correlation rejection, not a manufacturer safety limit."""
    def __init__(self, report):
        self.report=report
        super().__init__("Temperature is outside the positive-capacity material domain: "
                         f"{report['material']}, species_index={report['species_index']}, "
                         f"site_index={report['site_index']}, T={report['temperature_k']:.9g} K; "
                         f"qmax={report['capacity_mol_kg']:.9g} mol/kg")


def update_trajectory_extrema(extrema, accepted, bed, species_names=None):
    """Aggregate accepted BDF-node T/P ranges by actual adsorbent layer.

    ``accepted`` is (beds, cells, state variables, accepted times). This is
    an observational screen, not a continuous-time enclosure or fit proof.
    """
    # Local import avoids core -> qualification -> core import cycle.
    from .core import R, BAR
    if species_names is None:
        species_names = [f"species_{i+1}" for i in range(bed.ns)]
    species_names = list(species_names)
    if len(species_names) != bed.ns:
        raise ValueError("species_names must contain one name per gas species")
    for index, material in enumerate(bed.materials):
        cells = [i for i, assigned in enumerate(bed.layer_indices) if assigned is material]
        if not cells:
            raise ValueError(f"No cells allocated to material index {index}")
        t = accepted[:, cells, -2, :]
        if bed.equation_of_state == "peng_robinson":
            pressure_samples = []
            partial_samples = []
            for bed_index in range(accepted.shape[0]):
                bed_pressures = []
                bed_partials = []
                for time_index in range(accepted.shape[-1]):
                    fields = bed.fields(accepted[bed_index, :, :, time_index].ravel())
                    bed_pressures.append(fields[4][cells])
                    bed_partials.append(fields[4][cells, None] * fields[5][cells])
                pressure_samples.append(np.stack(bed_pressures, axis=-1))
                bed_partial = np.stack(bed_partials, axis=-1)
                partial_samples.append(bed_partial)
            p = np.stack(pressure_samples, axis=0)
            partial = np.stack(partial_samples, axis=0)
        else:
            gas = np.maximum(accepted[:, cells, :bed.ns, :], 1e-30)
            p = gas.sum(axis=2) * R * t / bed.vgas[cells][None, :, None] / BAR
            partial = gas * R * t[:, :, None, :] / bed.vgas[cells][None, :, None, None] / BAR
        if not np.isfinite(t).all() or not np.isfinite(p).all():
            raise FloatingPointError("Non-finite trajectory domain sample")
        key = str(index)
        previous = extrema.get(key, {})
        partial_ranges = {}
        previous_partial = previous.get("partial_pressure_bar", {})
        for species_index, species_name in enumerate(species_names):
            species_previous = previous_partial.get(species_name, {})
            values = partial[:, :, species_index, :]
            partial_ranges[species_name] = {
                "min": min(float(values.min()), species_previous.get("min", math.inf)),
                "max": max(float(values.max()), species_previous.get("max", -math.inf)),
            }
        extrema[key] = {
            "material": material["name"],
            "temperature_min_k": min(float(t.min()), previous.get("temperature_min_k", math.inf)),
            "temperature_max_k": max(float(t.max()), previous.get("temperature_max_k", -math.inf)),
            "pressure_min_bar": min(float(p.min()), previous.get("pressure_min_bar", math.inf)),
            "pressure_max_bar": max(float(p.max()), previous.get("pressure_max_bar", -math.inf)),
            "partial_pressure_bar": partial_ranges,
            "accepted_cell_time_samples": int(previous.get("accepted_cell_time_samples", 0) + t.size),
        }


def trajectory_domain_screen(input_screen, extrema, run_status):
    """Compare observed layer ranges with exact-grade declared fit ranges."""
    rows = []
    for index, material in enumerate(input_screen["materials"]):
        observed = extrema.get(str(index), {})
        documented = bool(material["experimental_fit_domain_documented"])
        within = None
        pressure_checks = None
        if documented and observed.get("accepted_cell_time_samples", 0) > 0:
            domain = material["experimental_fit_domain"]
            try:
                tlo, thi = map(float, domain["temperature_k"])
                temperature_within = (
                    tlo <= observed["temperature_min_k"]
                    <= observed["temperature_max_k"] <= thi
                )
                component_domain = domain.get("component_partial_pressure_bar")
                if component_domain:
                    observed_partial = observed.get("partial_pressure_bar", {})
                    pressure_checks = {}
                    for species_name, declared_range in component_domain.items():
                        plo, phi = map(float, declared_range)
                        sampled = observed_partial.get(species_name)
                        component_within = bool(
                            sampled
                            and plo <= float(sampled["min"])
                            <= float(sampled["max"]) <= phi
                        )
                        pressure_checks[species_name] = {
                            "documented_partial_pressure_bar": [plo, phi],
                            "observed_partial_pressure_bar": sampled,
                            "within_documented_domain": component_within,
                        }
                    pressure_within = (
                        set(pressure_checks) == set(observed_partial)
                        and all(row["within_documented_domain"] for row in pressure_checks.values())
                    )
                else:
                    plo, phi = map(float, domain["pressure_bar"])
                    pressure_within = (
                        plo <= observed["pressure_min_bar"]
                        <= observed["pressure_max_bar"] <= phi
                    )
                within = bool(temperature_within and pressure_within)
            except (KeyError, TypeError, ValueError):
                within = False
        rows.append({"material_index": index, "material": material["material"],
                     "experimental_fit_domain_documented": documented,
                     "documented_temperature_k": (material["experimental_fit_domain"] or {}).get("temperature_k"),
                     "documented_pressure_bar": (material["experimental_fit_domain"] or {}).get("pressure_bar"),
                     "documented_component_partial_pressure_bar": (
                         material["experimental_fit_domain"] or {}
                     ).get("component_partial_pressure_bar"),
                     "observed": observed or None,
                     "component_partial_pressure_checks": pressure_checks,
                     "accepted_nodes_within_documented_domain": within})
    all_observed = all((row["observed"] or {}).get("accepted_cell_time_samples", 0) > 0 for row in rows)
    all_within = all_observed and all(row["accepted_nodes_within_documented_domain"] is True for row in rows)
    return {"schema": "dynamic-psa-trajectory-domain-screen/v1",
            "case_id": input_screen["case_id"], "run_status": run_status,
            "sampling_scope": "All accepted BDF nodes through completed slots; not a continuous-time mathematical enclosure",
            "layers": rows, "all_material_layers_observed": all_observed,
            "accepted_nodes_all_within_documented_domains": all_within,
            "source_material_validation_released": False,
            "design_eligible": False, "training_eligible": False,
            "interpretation": "An input feed/ambient check cannot substitute for the simulated bed trajectory. Unknown fit domains or out-of-range accepted nodes remain fail-closed; passing this sampled check alone never qualifies the material or process."}


def check_capacity(material, temperature, capacity):
    """Report the first rejected cell/species/site, preserving complex tests."""
    cap=np.real(capacity)
    bad=np.argwhere(cap<0.)
    if not len(bad):
        return
    index=tuple(bad[0])
    species,site=index[-2:]
    t=np.real(np.asarray(temperature)).ravel()
    cell=int(index[0]) if cap.ndim==3 else None
    a1=np.asarray(material["a1"])
    a2=np.asarray(material.get("a2",np.zeros_like(a1)))
    slope=float(a2[species,site])
    raise MaterialDomainError({"category":"constitutive_positive_capacity_rejection",
        "material":material.get("name","unnamed"),"local_material_cell_index":cell,
        "species_index":int(species),"site_index":int(site),
        "temperature_k":float(t[cell or 0]),"capacity_mol_kg":float(cap[index]),
        "mathematical_zero_capacity_temperature_k":float(a1[species,site]/slope) if slope else None,
        "interpretation":"Mathematical validity of qmax=a1-a2*T only; neither a measured fit-domain limit nor a manufacturer safety rating."})


def case_screen(case):
    """Audit declared inputs without guessing unknown property-fit domains."""
    # Local imports avoid a core -> qualification -> core import cycle.
    from .core import Bed,valve_molar_flow
    from .schedules import build_schedule
    bed=Bed(case["bed"],case["species"])
    schedule=build_schedule(case["schedule"])
    species=[s["name"] for s in case["species"]]
    material_rows=[]
    for material in bed.materials:
        if material.get("model") == "extended_langmuir_3":
            a1=np.asarray(material["ip1_kmol_kg"],float) * 1000.0
            a2=np.asarray(material["ip2_kmol_kg_k"],float) * 1000.0
            equilibrium_model = "extended_langmuir_3"
        elif material.get("model") == "extended_langmuir_freundlich":
            # The source LF form has qmax(T)=k1+k2*T.  Its positive-capacity
            # screen is kept separate from the dual-site Ahn correlation.
            a1=np.asarray(material["lf_k1_mol_kg"],float)
            a2=-np.asarray(material["lf_k2_mol_kg_k"],float)
            equilibrium_model = "extended_langmuir_freundlich"
        elif material.get("model") == "streb_sips_wilson_rast_candidate":
            # Sips q_inf=A1*exp(A2*(T/Tref-1)) remains positive for any finite
            # A1>0; the unknown source Lambda matrix is an *independent*
            # scientific gate, not a mathematical capacity limit.
            a1=np.ones((len(species),1),float)
            a2=np.zeros_like(a1)
            equilibrium_model = "streb_sips_wilson_rast_candidate"
        elif material.get("model") == "guderian_2021_sips_iast":
            # The source evaluates temperature-dependent pure Sips models
            # through IAST.  Positivity and domain enforcement live in the
            # source-local constitutive module, not in legacy a1-a2*T arrays.
            a1=np.ones((len(species),1),float)
            a2=np.zeros_like(a1)
            equilibrium_model = "guderian_2021_sips_iast"
        elif material.get("model") == "competitive_sips_source":
            # Schell's source form uses a temperature-dependent Sips capacity
            # and affinity, rather than the engine's dual-site qmax=a1-a2*T
            # law.  Keep its positivity screen and fit-domain identity in a
            # separate branch so a source Sips table can never be interpreted
            # as an Ahn/Wang dual-site table.
            a1=np.asarray(material["sips_a_mol_kg"],float)
            a2=np.zeros_like(a1)
            equilibrium_model = "competitive_sips_source"
        elif material.get("model") == "ko_extended_langmuir_source":
            # Ko (2021) stores qsat=qsa+qsb/T rather than qmax=a1-a2*T.
            # Use a neutral mathematical positivity screen here; the actual
            # source-form capacity check is performed by core.equilibrium.
            a1=np.ones((len(species),1),float)
            a2=np.zeros_like(a1)
            equilibrium_model = "ko_extended_langmuir_source"
        elif material.get("model") == "dantas_toth_source":
            # Dantas' independent Toth capacities are not the common
            # dual-site qmax=a1-a2*T law.  Use qm only for a mathematical
            # non-negative-capacity screen and keep the source fit domain
            # closed until digitized breakthrough data are bound.
            a1=np.asarray(material["dantas_qm_mol_kg"], float)
            a2=np.zeros_like(a1)
            equilibrium_model = "dantas_toth_source"
        elif material.get("model") == "tamof1_dual_site_langmuir_freundlich_303k":
            # TAMOF-1 is a source-local dual-site Langmuir--Freundlich
            # namespace.  Its published 303 K capacities are q_m, not the
            # legacy Ahn/Wang a1-a2*T coefficients.  Use q_m only for the
            # mathematical non-negative-capacity screen and retain the
            # source model identity in the audit.  In particular, do not
            # manufacture a temperature slope or a documented dynamic fit
            # domain from this equilibrium table.
            a1=np.asarray(material["tamof_qm_mol_kg"], float)
            a2=np.zeros_like(a1)
            equilibrium_model = "tamof1_dual_site_langmuir_freundlich_303k"
        else:
            a1=np.asarray(material["a1"],float)
            a2=np.asarray(material.get("a2",np.zeros_like(a1)),float)
            equilibrium_model = material.get("equilibrium_model", "dual_site_langmuir")
        limits=[]
        for i,j in np.argwhere(a2>0.):
            limits.append({"species":species[i],"site_index":int(j),"zero_capacity_temperature_k":float(a1[i,j]/a2[i,j])})
        validity=material.get("experimental_fit_domain",{}) or {}
        if not isinstance(validity,dict):
            validity={}
        # A DOI and a temperature interval alone are not enough: the material
        # grade, species coverage and equilibrium model must be bound to the
        # exact case before a source-domain gate can open.  This deliberately
        # rejects attractive but incomplete literature metadata.
        required_domain_keys=("temperature_k","source","locator",
                              "material_grade","components","equilibrium_model")
        declared_components={str(name) for name in validity.get("components",[]) }
        case_components=set(species)
        case_grade=material.get("material_grade")
        grade_match=(case_grade not in (None, "") and validity.get("material_grade")==case_grade)
        model_match=validity.get("equilibrium_model")==equilibrium_model
        total_pressure_domain = validity.get("pressure_bar")
        component_pressure_domain = validity.get("component_partial_pressure_bar")
        component_pressure_keys = (
            {str(name) for name in component_pressure_domain}
            if isinstance(component_pressure_domain, dict) else set()
        )
        pressure_domain_declared = bool(
            total_pressure_domain not in (None, "", [])
            or (
                isinstance(component_pressure_domain, dict)
                and component_pressure_domain
                and component_pressure_keys == case_components
            )
        )
        known=(all(validity.get(k) not in (None,"",[]) for k in required_domain_keys)
               and pressure_domain_declared and declared_components==case_components
               and grade_match and model_match)
        operating_temperatures = [
            float(case["feed"]["temperature_k"]),
            float(case["bed"].get("ambient_k", case["feed"]["temperature_k"])),
        ]
        operating_pressures = [
            float(case["feed"]["pressure_bar"]),
            float(case.get("high_pressure_bar", case["feed"]["pressure_bar"])),
            float(case.get("product_pressure_bar", case["feed"]["pressure_bar"])),
        ]
        range_check = None
        operating_partial_pressures = None
        if known:
            try:
                t_range = [float(value) for value in validity["temperature_k"]]
                p_range = (
                    [float(value) for value in total_pressure_domain]
                    if total_pressure_domain not in (None, "", []) else []
                )
            except (TypeError, ValueError):
                t_range, p_range = [], []
            temperature_check = (
                len(t_range) == 2
                and min(t_range) <= min(operating_temperatures)
                and max(operating_temperatures) <= max(t_range)
            )
            if component_pressure_domain:
                feed_y = np.asarray(case["feed"]["mole_fractions"], dtype=float)
                initial_y = np.asarray(
                    case.get("initial_mole_fractions", [1.0] + [0.0] * (len(species) - 1)),
                    dtype=float,
                )
                boundary_partials = np.vstack((
                    feed_y * float(case["feed"]["pressure_bar"]),
                    initial_y * float(case.get("initial_pressure_bar", case["feed"]["pressure_bar"])),
                ))
                operating_partial_pressures = {
                    name: [float(boundary_partials[:, i].min()), float(boundary_partials[:, i].max())]
                    for i, name in enumerate(species)
                }
                try:
                    pressure_check = all(
                        len(component_pressure_domain[name]) == 2
                        and float(component_pressure_domain[name][0])
                        <= operating_partial_pressures[name][0]
                        <= operating_partial_pressures[name][1]
                        <= float(component_pressure_domain[name][1])
                        for name in species
                    )
                except (KeyError, TypeError, ValueError):
                    pressure_check = False
            elif len(p_range) == 2:
                pressure_check = (
                    min(p_range) <= min(operating_pressures)
                    and max(operating_pressures) <= max(p_range)
                )
            else:
                pressure_check = False
            if len(t_range) != 2:
                range_check = False
            else:
                range_check = bool(temperature_check and pressure_check)
        material_rows.append({"material":material["name"],"mathematical_positive_capacity_limits":limits,
            "equilibrium_model":equilibrium_model,
            "experimental_fit_domain":validity or None,"experimental_fit_domain_documented":known,
            "experimental_fit_domain_required_fields":list(required_domain_keys),
            "experimental_fit_domain_component_match":declared_components==case_components,
            "case_material_grade":case_grade,
            "experimental_fit_domain_grade_match":grade_match,
            "experimental_fit_domain_equilibrium_model_match":model_match,
            "operating_temperatures_k":operating_temperatures,
            "operating_pressures_bar":operating_pressures,
            "operating_boundary_partial_pressures_bar":operating_partial_pressures,
            "operating_point_within_documented_domain":range_check,
            "caveat":"Positivity is necessary but does not establish experimental validity. A different grade or binder is a different material."})
    hydraulic={"mode":case.get("boundary_mode"),"plant_valve_rating_established":False}
    if case.get("boundary_mode")=="compressible_valve" and "area_cd_m2" in case["valves"]["feed"]:
        feed=case["feed"]
        area=case["valves"]["feed"]["area_cd_m2"]
        upper=float(valve_molar_flow(feed["pressure_bar"],case["product_pressure_bar"],feed["temperature_k"],
            np.asarray(feed["mole_fractions"]),bed.mw,bed.cp,area,case.get("valve_regularization_bar",1e-6)))
        required=feed["molar_flow_mol_s"]
        hydraulic.update({"effective_area_cd_m2":area,"per_active_bed_feed_cap_mol_s":required,
            "optimistic_feed_limit_during_product_export_mol_s":upper,
            "instantaneous_reference_cap_passes_necessary_screen":upper>=required,
            "necessary_effective_area_lower_bound_m2":area*required/upper if upper>0. else None,
            "scope":"Forward-flow product-export periods only, ignoring additional bed/pipe losses. Not a vendor Cv or proof of average plant-feed feasibility."})
    evidence=case.get("engineering_evidence",{})
    open_gates=[]
    thermal_closure=case["bed"].get("thermal_closure",{}) or {}
    thermal_geometry={"wall_coefficient_basis":thermal_closure.get("wall_coefficient_basis"),
        "reference_lab_diameter_m":thermal_closure.get("reference_lab_diameter_m"),
        "case_diameter_m":float(case["bed"]["diameter_m"]),
        "radial_temperature_model":thermal_closure.get("radial_temperature_model"),
        "case_wall_coefficient_w_m2_k":float(case["bed"]["h_bed_wall_w_m2_k"]),
        "geometry_consistent_k_over_radius_w_m2_k":float(case["bed"]["bed_k_w_m_k"]/(case["bed"]["diameter_m"]/2.0)),
        "geometry_consistency_only_not_source_validation":True}
    if thermal_closure.get("wall_coefficient_basis")=="k_over_radius_screening_heuristic":
        if not math.isclose(thermal_geometry["case_wall_coefficient_w_m2_k"],
                            thermal_geometry["geometry_consistent_k_over_radius_w_m2_k"],rel_tol=1e-9):
            open_gates.append("wall_coefficient_inconsistent_with_case_geometry")
        radial_evidence=evidence.get("plant_scale_radial_heat_transfer",{})
        if not isinstance(radial_evidence,dict) or not all(radial_evidence.get(k) for k in
                ("source","locator","applicability","validated_diameter_range_m","temperature_gradient_evidence")):
            open_gates.append("plant_scale_radial_heat_transfer_validation")
        else:
            diameter_range=radial_evidence["validated_diameter_range_m"]
            try:
                diameter_in_range=(len(diameter_range)==2 and
                    float(diameter_range[0])<=thermal_geometry["case_diameter_m"]<=float(diameter_range[1]))
            except (TypeError,ValueError):
                diameter_in_range=False
            if not diameter_in_range:
                open_gates.append("plant_scale_radial_heat_transfer_validation")
    if schedule.get("name") == "streb_two_bed_8operation_screen":
        # Port topology and timing are executable, but this adapter uses
        # caller-supplied bottle/valve values and does not implement the
        # article's pressure trajectory or source-fitted Wilson/RAST closure.
        # No amount of engineering paperwork may promote it to a source fit.
        open_gates.append("streb_lab_topology_only_not_source_equivalent")
    if case.get("boundary_mode") == "streb_pressure_equation_screen":
        # The conservative flow inversion follows Eq. (3) for a caller-
        # supplied asymptote, but step-transition pressures and the meaning
        # of published P_low remain unresolved.  This is not source replay.
        open_gates.append("streb_pressure_equation_boundary_diagnostic_not_source_equivalent")
    if any(m.get("model") == "streb_sips_wilson_rast_candidate" for m in bed.materials):
        open_gates.append("streb_wilson_lambda_source_parameters_unresolved")
    if schedule.get("name") == "guderian_two_bed_6step_source":
        # The product receiver, four global intervals and the source Eq. (3)
        # velocity/pressure/temperature-dependent axial dispersion are now
        # executable.  PR departure internal energy and an independently
        # documented initial/CSS state remain separate source-equivalence
        # gates and therefore keep this adapter fail closed.
        if not bed.pr_departure_energy:
            open_gates.append("guderian_pr_departure_energy_not_implemented")
        open_gates.append("guderian_initial_state_and_published_css_norm_unavailable")
    if not schedule["timing_inputs_unambiguous"]:
        open_gates.append("schedule_timing_inputs_ignored_by_plugin")
    if not schedule["cycle_time_consistent"]:
        open_gates.append("declared_cycle_time_conflicts_with_executed_schedule")
    if not all(row["experimental_fit_domain_documented"] for row in material_rows):
        open_gates.append("experimental_isotherm_temperature_pressure_and_material_grade_domain")
    if any(row["operating_point_within_documented_domain"] is False for row in material_rows):
        open_gates.append("declared_operating_point_outside_documented_material_domain")
    for key in ("thermal_properties","layer_density_and_heat_capacity","valves_and_pressure_loss_budget",
                "vessel_mechanical_design","drying_guard_and_regeneration","upstream_average_feed_and_delivery_pressure"):
        record=evidence.get(key,{})
        if not isinstance(record,dict) or not all(record.get(k) for k in ("source","locator","applicability")):
            open_gates.append(key)
    return {"schema":"dynamic-psa-engineering-input-screen/v1","case_id":case["case_id"],
        "number_of_beds":schedule["number_of_beds"],"number_of_cycle_steps":schedule["number_of_cycle_steps"],
        "schedule_timing_inputs_used":schedule["used_timing_inputs"],
        "schedule_timing_inputs_ignored":schedule["ignored_timing_inputs"],
        "total_adsorbent_mass_kg":float(bed.mass.sum()*schedule["number_of_beds"]),
        "total_bed_volume_m3":float(bed.cell_volume.sum()*schedule["number_of_beds"]),
        "materials":material_rows,"hydraulic_necessary_screen":hydraulic,
        "thermal_geometry_screen":thermal_geometry,
        "equation_of_state":bed.equation_of_state,"open_engineering_evidence_gates":open_gates,
        "engineering_evidence_complete":not open_gates,"training_eligible":False,
        "interpretation":"Input evidence screen only. CSS, conservation, external validation, independent resolution and average-flow matching remain separate mandatory gates."}
