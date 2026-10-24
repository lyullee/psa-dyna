"""Traceable source and design-case factories. No expensive imports or I/O."""
from copy import deepcopy
import hashlib
import json
from .core import R, BAR


def wang_case(schedule="wang_type_i_12step", n_cells=4):
    """Published RB3 apparatus; closure parameters remain explicitly declared."""
    material = {
        "name":"Wang2023_RB3_C1", "end_fraction":1.0,
        "material_grade":"Norit RB3 (Sigma-Aldrich)",
        "a1":[[2.97,0],[4.56,2.30],[2.26,0],[0.75,7.69]],
        "b0_per_bar":[[2.37e-4,0],[3.03e-5,6.77e-4],[2.67e-4,0],[4.19e-5,3.50e-5]],
        "binding_energy_j_mol":[[9540,0],[18570,16860],[16360,0],[29990,22380]],
        "ldf_per_s":[2,3,3,1],
        # Wang et al. Figure 5 reports pure-component isotherms on this exact
        # commercial grade.  H2/CH4 were measured at 273.15, 303.15 and
        # 333.15 K; N2/CO2 additionally include 293.15 and 313.15 K.  The
        # pressure axes are component-specific, so qualification must compare
        # simulated partial pressures rather than total bed pressure.
        "experimental_fit_domain":{
            "temperature_k":[273.15,333.15],
            "component_partial_pressure_bar":{
                "H2":[0.0,7.0], "CH4":[0.0,12.0],
                "N2":[0.0,1.2], "CO2":[0.0,1.2],
            },
            "measured_temperature_k":{
                "H2":[273.15,303.15,333.15],
                "CH4":[273.15,303.15,333.15],
                "N2":[273.15,293.15,313.15,333.15],
                "CO2":[273.15,293.15,313.15,333.15],
            },
            "source":"10.1016/j.cej.2023.144920",
            "locator":"Article Sections 2.1 and 3.1, Figure 5 and Table 3",
            "material_grade":"Norit RB3 (Sigma-Aldrich)",
            "components":["H2","CH4","N2","CO2"],
            "equilibrium_model":"dual_site_langmuir",
            "pressure_basis":"pure_component_partial_pressure",
        },
    }
    case = {
        "schema":"dynamic-psa-case/v1", "case_id":f"wang_{schedule}_{n_cells}cell",
        "species":[{"name":name,"mw_kg_mol":mw,"cp_j_mol_k":cp} for name,mw,cp in
                   zip(["H2","CH4","N2","CO2"],[.00201588,.0160425,.0280134,.0440095],[28.84,35.69,29.12,37.13])],
        "bed":{"n_cells":n_cells,"length_m":1.07,"diameter_m":.02118,"void_fraction":.337,
               "pellet_void_fraction":.65,"bulk_density_kg_m3":470.,"particle_diameter_m":.00276,
               "solid_cp_j_kg_k":1046.,"bed_k_w_m_k":.24,"wall_thickness_m":.00211,
               "wall_density_kg_m3":7850.,"wall_cp_j_kg_k":460.,"wall_k_w_m_k":10.9,
               "h_bed_wall_w_m2_k":.24/(.02118/2),"h_ambient_w_m2_k":.006,
               "ambient_k":293.95,"kinetics":"partial_pressure_ldf","energy_mode":"source_sensible",
               "materials":[material]},
        "feed":{"mole_fractions":[.58,.29,.07,.06],"temperature_k":293.95,
                "pressure_bar":5.,"molar_flow_mol_s":.00223},
        "high_pressure_bar":5.,"low_pressure_bar":.1,"product_pressure_bar":5.,
        "initial_pressure_bar":.1,"purge_to_feed_ratio":.1,"purge_gain_mol_s_bar":100.,
        "boundary_mode":"wang_controller",
        "valve_regularization_bar":1e-6,
        "valves":{"feed":{"taper_band_bar":.35},"product":{"gain_mol_s_bar":.02},
                  "tail":{"gain_mol_s_bar":.01},"equalization":{"gain_mol_s_sqrt_bar":.003}},
        # Populate only timing keys consumed by the selected plugin.  Older
        # archived Wang cases also carried pressurization_s=30 even though
        # the 12-slot plugin never applied it; those archives remain intact.
        "schedule":{"name":schedule},
        "solver":{"rtol":2e-7,"max_step_s":8.,"maximum_cycles":150,"minimum_cycles":3,
                  "css_tolerance":1e-4,"css_consecutive_cycles":3,
                  "mass_balance_tolerance":1e-5,"energy_balance_tolerance":1e-5},
        "provenance":{"doi":"10.1016/j.cej.2023.144920","material_basis":"Table3/S3",
                      "geometry_basis":"Tables1/2","boundary_closures":["valve gains","bed-wall h=k/r","equal-slot timing"],
                      "energy_mode":"Published solid-dominated thermal closure; gas internal-energy terms omitted explicitly."},
        "claim":"Source-apparatus reproduction, not validation of S-015 design."
    }
    if schedule in {"breakthrough", "closed"}:
        case["schedule"]["duration_s"] = 1600.0
    elif schedule == "skarstrom_4step":
        case["schedule"].update({"pressurization_s":30., "adsorption_s":120.})
    else:
        case["schedule"]["slot_s"] = 120.
    if schedule == "breakthrough":
        case["initial_pressure_bar"] = 5.
        case["solver"].update({"record_trace":True,"trace_interval_s":2.})
    else:
        # The local source package contains a 1600-s continuous breakthrough
        # duration.  It is retained as explicit provenance, not as a claimed
        # cyclic period; the published Type-I valve timing is inferred from
        # the topology/pressure trace until a complete valve log is found.
        case["schedule"]["reported_breakthrough_duration_s"] = 1600.0
        # Table 4 annotates Run XII's AD-flow cell as 2.6 SLPM (40.5 degC).
        # We apply 40.5 C at the feed boundary, but the table does not state
        # where that temperature was measured; retain the assumption and do
        # not treat it as a source-verified bed-inlet sensor reading.
        case["feed"].update({"pressure_bar":5.17,"molar_flow_mol_s":2.6/3*.00223,"temperature_k":313.65})
        case["provenance"]["cycle_feed_temperature_basis"]="Main article Table 4 Run XII AD-flow cell: 2.6 SLPM (40.5 degC); modeled as bed-inlet gas T, but measurement location is not stated."
        case["provenance"]["cycle_feed_temperature_location_verified"] = False
        case["provenance"]["cycle_feed_basis"]="adsorption_active_bed"
        case.update({"high_pressure_bar":5.17,"product_pressure_bar":5.17,
                     "low_pressure_bar":.11,"initial_pressure_bar":.11,"purge_to_feed_ratio":.0077})
    return case


def dantas_2011_co2_n2_fixed_bed_case(run=5, n_cells=16):
    """Build a source-local Dantas fixed-bed diagnostic case.

    Runs 1--4 are CO2/He and runs 5--8 are CO2/N2.  The common engine keeps
    helium as an explicit inert component rather than silently replacing it
    with nitrogen.  Runs 9--12 are intentionally rejected here because the
    published non-adiabatic external wall closure is not complete enough for
    a source-equivalent replay.  The result is a single-bed breakthrough
    case, never a cyclic VPSA or design case.
    """
    import math
    from .dantas_2011_co2_n2_adapter import DantasCO2N2Zeolite13X, DantasCase, SOURCE_ID

    run = int(run)
    if run not in range(1, 9):
        raise ValueError("Dantas fixed-bed engine cases are limited to source runs 1--8; runs 9--12 remain closed")
    if int(n_cells) < 2:
        raise ValueError("n_cells must be at least 2")
    source = DantasCO2N2Zeolite13X()
    metadata = DantasCase(run).manifest()
    components = ("CO2", "N2", "He")
    material = source.build_engine_material(components)
    temperature = float(metadata["temperature_K"])
    pressure = float(metadata["pressure_bar"])
    flow_m3_s = float(metadata["flow_m3_s"])
    molar_flow = flow_m3_s * pressure * BAR / (R * temperature)
    feed = [0.20, 0.0, 0.80] if run <= 4 else [0.20, 0.80, 0.0]
    ldf = metadata["ldf_overall_s_inv"]
    # Loading-LDF kinetics are represented by one coefficient per gas
    # component.  Keep this as a flat vector: the common Bed validator
    # deliberately rejects nested species-by-site arrays for this mode.
    material["ldf_per_s"] = [float(ldf.get(name, 0.0)) for name in components]
    length = 0.171
    diameter = 0.022
    volume = math.pi * diameter**2 / 4.0 * length
    bed_mass = 0.0455
    case = {
        "schema": "dynamic-psa-case/v1",
        "case_id": f"dantas_2011_run_{run}_co2_fixed_bed_{int(n_cells)}cell",
        "source_id": SOURCE_ID,
        "adsorbent_id": "zeolite_13X",
        "reporting_species": "CO2",
        "species": [
            {"name": "CO2", "mw_kg_mol": 0.0440095, "cp_j_mol_k": 37.13},
            {"name": "N2", "mw_kg_mol": 0.0280134, "cp_j_mol_k": 29.12},
            {"name": "He", "mw_kg_mol": 0.0040026, "cp_j_mol_k": 20.79},
        ],
        "bed": {
            "n_cells": int(n_cells), "length_m": length, "diameter_m": diameter,
            "void_fraction": 0.43, "pellet_void_fraction": 0.37,
            "bulk_density_kg_m3": bed_mass / volume, "particle_diameter_m": 0.0029,
            "solid_cp_j_kg_k": 920.0, "bed_k_w_m_k": 0.20,
            # Dantas Table 4: D_L = 0.08 x 10^-4 m2/s for runs 1--4 and
            # 0.15 x 10^-4 m2/s for runs 5--8. This is source-local axial
            # dispersion, not a Wang/biogas transfer.
            "axial_dispersion_m2_s": (0.08e-4 if run <= 4 else 0.15e-4),
            "wall_thickness_m": 0.0015, "wall_density_kg_m3": 7830.0,
            "wall_cp_j_kg_k": 440.0, "wall_k_w_m_k": 13.4,
            "h_bed_wall_w_m2_k": 0.0, "h_ambient_w_m2_k": 0.0,
            "ambient_k": temperature, "kinetics": "loading_ldf",
            "energy_mode": "internal_energy", "materials": [material],
        },
        "feed": {"mole_fractions": feed, "temperature_k": temperature,
                 "pressure_bar": pressure, "molar_flow_mol_s": molar_flow},
        # The source experiment is approximately isobaric.  The engine still
        # requires a strictly ordered high/low pair for its valve contract;
        # the 0.01-bar low value is a numerical floor, while the isobaric
        # boundary controller holds the actual bed pressure at the source P.
        "high_pressure_bar": pressure, "low_pressure_bar": 0.01,
        "product_pressure_bar": pressure, "initial_pressure_bar": pressure,
        "initial_mole_fractions": [0.0, 0.0, 1.0],
        "purge_to_feed_ratio": 0.0, "purge_gain_mol_s_bar": 0.0,
        "boundary_mode": "isobaric_breakthrough", "pressure_control_gain_mol_s_bar": 1.0e-3,
        "valve_regularization_bar": 1.0e-6,
        "valves": {"feed": {"taper_band_bar": 0.05}, "product": {"gain_mol_s_bar": 0.0},
                    "tail": {"gain_mol_s_bar": 0.0}, "equalization": {"gain_mol_s_sqrt_bar": 0.0}},
        "schedule": {"name": "breakthrough", "duration_s": 36000.0 if run <= 4 else 18000.0},
        "solver": {"rtol": 2e-7, "max_step_s": 8.0, "maximum_cycles": 1,
                    "minimum_cycles": 1, "css_tolerance": 1e-4,
                    "css_consecutive_cycles": 1, "mass_balance_tolerance": 1e-5,
                    "energy_balance_tolerance": 1e-5, "record_trace": True,
                    "trace_interval_s": 5.0},
        "provenance": {
            "source_doi": "10.1590/S0104-66322011000300018",
            "source_run": run, "source_case": metadata,
            "flow_basis": "published volumetric flow converted with ideal gas at source P/T",
            "initial_state_basis": "He-regenerated bed represented explicitly as the initial carrier gas",
            "wall_boundary_basis": "adiabatic diagnostic closure for runs 1--8; wall exchange coefficients set to zero",
            "parameter_isolation": SOURCE_ID,
        },
        "claim": "Dantas source-local CO2 fixed-bed diagnostic; figure-digitized experimental comparison pending; not cyclic CSS, DWSIM, TEA/LCA or design release.",
    }
    payload = {"source_id": SOURCE_ID, "run": run, "species": case["species"],
               "bed": case["bed"], "feed": case["feed"], "boundary_mode": case["boundary_mode"]}
    case["parameter_hash"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return case


def guderian_2021_n2_reference_case(experiment_id=13, n_cells=8):
    """Source-isolated CT-350 twin-bed N2-PSA aggregate-validation case.

    IDs 13--15 are the public workbook rows for which feed pressure, product
    pressure, all cycle timings and all three mean flow rates are populated.
    The factory deliberately does not fill blank pressures in other rows by
    relabelling product pressure.  The two physical beds execute six named
    operations over four simultaneous global time intervals.

    The current adapter includes the dynamic 12-L product receiver, published
    Sips/IAST/QDF laws, Peng--Robinson inventory, fixed flow controllers and
    source linear resistances.  Variable axial dispersion and PR departure
    energy remain explicit qualification gates, so this factory is not yet a
    source-equivalent validation result.
    """
    rows = {
        13: {"oxygen_ppm": 1000.0, "feed_pressure_bar": 8.45,
             "product_pressure_bar": 8.0, "feed_nm3_h": 1.335688016,
             "product_nm3_h": 0.4422808, "purge_nm3_h": 0.04498,
             "productivity_nm3_h_m3": 111.8, "air_demand": 3.02},
        14: {"oxygen_ppm": 100.0, "feed_pressure_bar": 8.45,
             "product_pressure_bar": 8.0, "feed_nm3_h": 1.1449653,
             "product_nm3_h": 0.2694036, "purge_nm3_h": 0.04498,
             "productivity_nm3_h_m3": 68.1, "air_demand": 4.25},
        15: {"oxygen_ppm": 10.0, "feed_pressure_bar": 8.45,
             "product_pressure_bar": 8.0, "feed_nm3_h": 1.00209436,
             "product_nm3_h": 0.1372732, "purge_nm3_h": 0.04498,
             "productivity_nm3_h_m3": 34.7, "air_demand": 7.3},
    }
    if experiment_id not in rows:
        raise ValueError("Guderian reference adapter accepts only source-complete experimental IDs 13, 14 or 15")
    if isinstance(n_cells, bool) or int(n_cells) != n_cells or n_cells < 2:
        raise ValueError("n_cells must be an integer of at least two")
    source = rows[experiment_id]
    temperature_k = 293.15

    def normal_molar_flow(value_nm3_h):
        return float(value_nm3_h) / 3600.0 * 101325.0 / (R * 273.15)

    def linear_series_gain_mol_s_bar(*coefficients_kmol_s_bar):
        # Linear pressure losses in series add as inverse conductances.
        coefficient = 1.0 / sum(1.0 / float(value) for value in coefficients_kmol_s_bar)
        return coefficient * 1000.0

    armature_gain = linear_series_gain_mol_s_bar(1.1e-4, 9.5e-4)
    # Article Table 4 identifies V1 as a linear control valve (the feed-flow
    # meter), not an imposed mass-flow controller.  The feed path also contains
    # V2/V3, the supporting plate and the microfilter.  Treating the workbook's
    # cycle-average feed flow as an instantaneous source kept adding gas above
    # the 8.45-bar supply boundary during repeated cycles.  The correct source
    # contract is the series linear conductance below; the workbook flow remains
    # an independent cycle-average validation target.
    feed_path_gain = linear_series_gain_mol_s_bar(9.0e-4, 5.0e-1, 1.1e-4, 9.5e-4)
    equalization_gain = linear_series_gain_mol_s_bar(
        2.5e-5, 1.1e-4, 9.5e-4, 1.1e-4, 9.5e-4
    )
    oxygen_product = source["oxygen_ppm"] * 1.0e-6
    material = {
        "name": "Shirasagi_MSC_CT-350_Guderian2021",
        "material_grade": "Shirasagi MSC CT-350",
        "model": "guderian_2021_sips_iast",
        "kinetic_model": "guderian_qdf_darken",
        "species_order": ["O2", "N2"],
        "parameter_provenance": "Marcinek et al. 2021, doi:10.1007/s10450-021-00320-0, Tables 1--3 and Eqs. 4--12",
        "end_fraction": 1.0,
        "experimental_fit_domain": {
            "temperature_k": [293.15, 333.15],
            "pressure_bar": [0.0, 20.05],
            "source": "Marcinek et al. 2021 and its CC-BY supplement",
            "locator": "Article Tables 1--3; supplement pure-isotherm Tables 18--19",
            "material_grade": "Shirasagi MSC CT-350",
            "components": ["O2", "N2"],
            "equilibrium_model": "guderian_2021_sips_iast",
        },
    }
    species = [
        {"name": "O2", "mw_kg_mol": 0.0319988, "cp_j_mol_k": 29.376,
         "critical_temperature_k": 154.581, "critical_pressure_pa": 5.043e6,
         "acentric_factor": 0.0222},
        {"name": "N2", "mw_kg_mol": 0.0280134, "cp_j_mol_k": 29.124,
         "critical_temperature_k": 126.192, "critical_pressure_pa": 3.3958e6,
         "acentric_factor": 0.0372},
    ]
    # Source left/right beds differ by only 2 mm and 8.47 g.  The current
    # common Network shares one Bed object; use the arithmetic mean geometry
    # and preserve the exact source values in provenance rather than claiming
    # heterogeneous-bed replay.
    length_m = (0.582 + 0.580) / 2.0
    packed_density_kg_m3 = (712.005 + 710.247) / 2.0
    case = {
        "case_id": f"Guderian2021_CT350_exp{experiment_id}_{int(n_cells)}cell",
        "species": species,
        "reporting_species": "N2",
        "bed": {
            "n_cells": int(n_cells), "length_m": length_m, "diameter_m": 0.066,
            "void_fraction": 0.404, "pellet_void_fraction": 0.234,
            # Eq. 13 uses 2*r_p/shape_factor as the Ergun hydraulic size.
            "particle_diameter_m": 2.0 * 0.830e-3 / 1.959,
            "bulk_density_kg_m3": packed_density_kg_m3,
            "solid_cp_j_kg_k": 880.0, "bed_k_w_m_k": 0.675,
            "wall_model": "fixed_ambient", "wall_thickness_m": 0.002,
            "h_bed_wall_w_m2_k": 50.0, "h_ambient_w_m2_k": 0.0,
            "ambient_k": temperature_k, "viscosity_pa_s": 1.81e-5,
            "kinetics": "loading_ldf", "energy_mode": "internal_energy",
            "equation_of_state": "peng_robinson",
            "pr_departure_energy": True,
            "pr_binary_interaction": [[0.0, 0.0], [0.0, 0.0]],
            "axial_dispersion_m2_s": 0.0,
            "axial_dispersion_model": "guderian_2021_eq3",
            "dispersion_particle_radius_m": 0.830e-3,
            "thermal_closure": {
                "wall_coefficient_basis": "published_fixed_ambient_gas_wall_coefficient",
                "reference_lab_diameter_m": 0.066,
                "radial_temperature_model": "one_dimensional_axial_fixed_ambient_wall",
            },
            "materials": [material],
        },
        "feed": {"mole_fractions": [0.209, 0.791],
                 "temperature_k": temperature_k,
                 "pressure_bar": source["feed_pressure_bar"],
                 "molar_flow_mol_s": normal_molar_flow(source["feed_nm3_h"])},
        "high_pressure_bar": source["feed_pressure_bar"],
        "low_pressure_bar": 1.01325,
        "product_pressure_bar": source["product_pressure_bar"],
        "initial_pressure_bar": 1.01325,
        "initial_mole_fractions": [0.209, 0.791],
        "initial_bed_pressures_bar": [source["product_pressure_bar"], 1.01325],
        "initial_bed_mole_fractions": [
            [oxygen_product, 1.0-oxygen_product], [0.209, 0.791]
        ],
        "product_receiver": {
            "volume_m3": 0.012,
            "initial_pressure_bar": source["product_pressure_bar"],
            "initial_temperature_k": temperature_k,
            "initial_mole_fractions": [oxygen_product, 1.0-oxygen_product],
            "thermal_mode": "adiabatic_sensible",
            "source_basis": "Marcinek dissertation reference receiver volume 12 L",
        },
        "boundary_mode": "guderian_receiver_cycle",
        "source_boundary": {
            "ambient_tail_pressure_bar": 1.01325,
            "purge_molar_flow_mol_s": normal_molar_flow(source["purge_nm3_h"]),
            "product_molar_flow_mol_s": normal_molar_flow(source["product_nm3_h"]),
            "observed_cycle_mean_feed_mol_s": normal_molar_flow(source["feed_nm3_h"]),
            "feed_boundary_mode": "linear_pressure_driven_from_8.45_bar_supply",
            "parameter_provenance": (
                "Public workbook IDs 13--15 plus Marcinek et al. 2021 Table 4; "
                "normal flow converted at 273.15 K and 1.01325 bar"
            ),
        },
        "valves": {
            "feed": {"gain_mol_s_bar": feed_path_gain},
            "product": {"gain_mol_s_bar": armature_gain},
            "tail": {"gain_mol_s_bar": armature_gain},
            "receiver_connection": {"gain_mol_s_bar": armature_gain},
            "equalization": {"gain_mol_s_bar": equalization_gain},
        },
        "valve_regularization_bar": 1.0e-8,
        "purge_to_feed_ratio": 0.0,
        "purge_gain_mol_s_bar": 0.0,
        "schedule": {"name": "guderian_two_bed_6step_source",
                     "half_cycle_s": 60.0, "equalization_s": 1.0,
                     "declared_cycle_time_s": 120.0},
        "solver": {"rtol": 2.0e-6, "max_step_s": 1.0,
                   "maximum_cycles": 80, "minimum_cycles": 5,
                   "css_tolerance": 2.0e-4, "css_consecutive_cycles": 3,
                   "mass_balance_tolerance": 2.0e-5,
                   "energy_balance_tolerance": 2.0e-5,
                   "record_trace": True, "trace_interval_s": 0.5,
                   "record_boundary_diagnostics": True,
                   "ledger_integration": "quadrature",
                   "ledger_epsabs": 1.0e-10, "ledger_epsrel": 1.0e-8,
                   "jacobian_method": "bounded_fd", "jacobian_relative_step": 1.0e-7,
                   "maximum_step_wall_s": 900.0},
        "engineering_evidence": {
            "thermal_properties": {"source": "Marcinek et al. 2021",
                "locator": "Table 3 and premises (5)--(6)",
                "applicability": "CT-350 reference apparatus only"},
            "layer_density_and_heat_capacity": {"source": "Marcinek et al. 2021 supplement",
                "locator": "Table 17 and article Table 3",
                "applicability": "mean of measured left/right packed beds"},
            "valves_and_pressure_loss_budget": {"source": "Marcinek et al. 2021",
                "locator": "Eq. 14 and Table 4",
                "applicability": "linear series conductance; ideal ball valves omitted"},
            "drying_guard_and_regeneration": {"source": "Marcinek et al. 2021 supplement",
                "locator": "dry-air apparatus description",
                "applicability": "source experiment only"},
            "upstream_average_feed_and_delivery_pressure": {"source": "Mendeley Data d2vx7w6x74.1",
                "locator": f"Experimental row ID {experiment_id}",
                "applicability": "cycle-mean source target"},
        },
        "provenance": {
            "source_doi": "10.1007/s10450-021-00320-0",
            "workbook_doi": "10.17632/d2vx7w6x74.1",
            "workbook_sha256": "2c9767b2423800e783dda01a9d319df0825233d70e10c311b2f860d26d713bc5",
            "experimental_row": {"id": experiment_id, **source},
            "exact_source_beds": {
                "left": {"length_m": 0.582, "volume_l": 1.990, "mass_kg": 1.41689},
                "right": {"length_m": 0.580, "volume_l": 1.983, "mass_kg": 1.40842},
            },
            "common_bed_approximation": "arithmetic mean geometry/density",
            "parameter_isolation": "O2/N2/CT-350 only; no H2 or CO2/CH4 parameters",
            "open_source_equivalence_gates": [
                "exact heterogeneous left/right bed geometry",
                "published initial state and numeric CSS norm",
            ],
        },
        "claim": (
            "Guderian/Marcinek source-local O2/N2 aggregate-validation candidate; "
            "not transient source equivalence and not transferable to H2 or biogas PSA."
        ),
    }
    payload = {
        "source_doi": case["provenance"]["source_doi"],
        "workbook_sha256": case["provenance"]["workbook_sha256"],
        "experiment_id": experiment_id, "species": species,
        "bed": case["bed"], "feed": case["feed"],
        "source_boundary": case["source_boundary"],
        "product_receiver": case["product_receiver"],
        "valves": case["valves"],
        "schedule": case["schedule"],
    }
    case["parameter_hash"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return case


def s015_case(schedule="skarstrom_4step", n_cells=8):
    """Dry S-015 composition; AC/5A design parameters are literature based.

    The supplied flow is the per-active-bed MFC cap. The simulation reports
    average external feed, so fluctuating simultaneous adsorption counts are
    not silently equated to the plant's continuous upstream delivery rate.
    """
    case = wang_case(schedule, n_cells)
    case["case_id"] = f"S015_dry_{schedule}_{n_cells}cell"
    order = ["H2","CO2","CH4","CO","N2"]
    case["species"] = [{"name":name,"mw_kg_mol":mw,"cp_j_mol_k":cp} for name,mw,cp in zip(
        order,[.00201588,.0440095,.0160425,.0280101,.0280134],[28.84,37.13,35.69,29.14,29.12])]
    def material(name, end, a1, a2, b0, b2, ldf, heat_cal_mol):
        return {"name":name,"end_fraction":end,"a1":[[v] for v in a1],"a2":[[v] for v in a2],
                "b0_per_bar":[[v] for v in b0],"affinity_energy_j_mol":[[v*R] for v in b2],
                "binding_energy_j_mol":[[v*4.184] for v in heat_cal_mol],"ldf_per_s":ldf}
    ac = material("Ahn_corrected_AC",.65,[16.943,28.7973,23.86,33.85,1.6441],
                  [.021,.07,.05621,.09072,.00073],[6.166e-5,.0099,.003432,.0002281,.0538],
                  [1229,1030,1159,1751,326],[.7,.0355,.195,.15,.261],[2880,5240,4290,4300,1660])
    zeo = material("Ahn_corrected_5A",1.,[4.314,10.03,5.833,11.8454,4.8133],
                   [.0106,.01858,.01192,.0313,.00668],[.002482,1.5575,.00059716,.0199,.000562],
                   [458,207,1731,763,1531],[.7,.0135,.147,.063,.099],[2800,9330,5400,5300,5470])
    # Wang's k/r screening coefficient belongs to a 21.18 mm lab tube.
    # Recompute the same *heuristic* for the S-015 diameter; this is not a
    # validated radial heat-transfer correlation for a 1.5 m packed bed.
    s015_diameter_m = 1.5
    s015_wall_screen_w_m2_k = case["bed"]["bed_k_w_m_k"] / (s015_diameter_m / 2.)
    case["bed"].update({"length_m":1.5,"diameter_m":s015_diameter_m,"void_fraction":.433,
                         "bulk_density_kg_m3":522*(1-.433),"pellet_void_fraction":0.,
                         "h_bed_wall_w_m2_k":s015_wall_screen_w_m2_k,
                         "thermal_closure":{
                             "wall_coefficient_basis":"k_over_radius_screening_heuristic",
                             "reference_lab_diameter_m":.02118,
                             "radial_temperature_model":"lumped_1d_axial",
                         },
                         "ambient_k":310.264937331605,"kinetics":"loading_ldf",
                         "energy_mode":"internal_energy","materials":[ac,zeo]})
    raw = [.7166470477,.2705743311,.0035538313,.0052723794,0.]
    case["feed"] = {"mole_fractions":[v/sum(raw) for v in raw],"temperature_k":310.264937331605,
                    "pressure_bar":15.6533970226832,"molar_flow_mol_s":25.1822807772193}
    case.update({"boundary_mode":"compressible_valve","high_pressure_bar":15.6533970226832,
                 "low_pressure_bar":1.1,"product_pressure_bar":15.,"initial_pressure_bar":1.1,
                 "purge_to_feed_ratio":.05,"purge_gain_mol_s_bar":10.,"valve_regularization_bar":1e-4})
    case["valves"]={"feed":{"area_cd_m2":2e-4},"product":{"area_cd_m2":2e-4},
                    "tail":{"area_cd_m2":4e-4},"equalization":{"area_cd_m2":2e-4}}
    if schedule == "skarstrom_4step":
        case["schedule"]={"name":schedule,"pressurization_s":30.,"adsorption_s":60.}
    else:
        case["schedule"]={"name":schedule,"slot_s":30.}
    # Plant-scale Ergun pressure differences are tiny relative to absolute
    # pressure. A smaller bounded Jacobian probe is a numerical setting,
    # not a change to the material data or governing equations.
    case["solver"]["jacobian_relative_step"]=1e-10
    case["provenance"]={"feed":"selected-feed forward DWSIM S-015, 2026-10-14",
                        "feed_water_mole_fraction_removed":.0039524105,
                        "material":"Ahn et al. research article DOI 10.1016/j.ces.2011.09.053 with its corrigendum DOI 10.1016/j.ces.2013.06.043; the corrigendum corrects the original Table 2 heat unit to cal/mol. Independent heat is converted with 4.184 J/cal; Affinity retains exp(k4/T).",
                        "engineering_assumptions":["bed diameter/length","valve areas","layer split","wall coefficient k/r recalculated for 1.5 m diameter; geometry-consistent heuristic only, not validated radial heat transfer","other thermal properties inherited as screening values"],
                        "water_guard":"Dry-feed design boundary. Water adsorption and guard regeneration/sizing remain separate."}
    case["claim"]="Literature-parameterized S-015 design screen; not a source-matched experiment or rated vessel."
    return deepcopy(case)


def bham_2019_co2_n2_breakthrough_case(n_cells=8, flow_basis="main_text_120_Nml_min"):
    """Independent CO2/N2 fixed-bed breakthrough adapter for Birmingham data.

    This is deliberately not a Wang/H2 or biogas case.  The source workbook
    reports a 40/60 mol% CO2/N2 feed at 25 bar and 298.15 K through an
    AC-MEA-MDEA bed.  Source equilibrium coefficients are converted from
    Pa^-1 to the engine's bar^-1 basis at the reference temperature.  The
    source reports one lumped CO2 LDF coefficient; the same coefficient is
    applied to the carrier N2 only as an explicitly flagged adapter
    assumption, so the run cannot be mistaken for a source-fitted validation
    until the N2 kinetic convention is independently resolved.
    """
    tref = 298.15
    e = [24800.0, 8400.0]
    b_ref_pa = [3.0e-6, 2.0e-10]
    b0_bar = [b * 1.0e5 * __import__("math").exp(-en / (R * tref)) for b, en in zip(b_ref_pa, e)]
    material = {
        "name": "AC-MEA-MDEA_Birmingham2019",
        "end_fraction": 1.0,
        "a1": [[9.2], [0.14]],
        "a2": [[0.0], [0.0]],
        "b0_per_bar": [[b0_bar[0]], [0.0]],
        "binding_energy_j_mol": [[e[0]], [e[1]]],
        "ldf_per_s": [0.046, 0.0],
        "parameter_provenance": "Birmingham 2019 Fuel 254 117492; machine-readable UBIRe eData record 344; CO2 K reported 0.046 s^-1. Article cites 0.14 mmol/g N2 capacity from a separate study but does not provide a matched N2 fit; N2 is therefore inert in this diagnostic adapter rather than assigned a fabricated kinetic coefficient.",
    }
    # The paper reports a measured 10 g sorbent charge.  The engine's
    # ``bulk_density_kg_m3`` multiplies the solid fraction of bed volume, so
    # convert that measured charge into an equivalent particle-density input
    # rather than silently trusting the table's separately printed 262 kg/m3
    # label.  The original value remains recorded in provenance.
    bed_volume = __import__("math").pi * 0.025**2 / 4.0 * 0.069
    # ``Bed.mass`` is defined per total packed-bed volume, so the measured
    # 10 g charge maps directly to 10 g / V_bed; the void fraction is already
    # represented separately in gas and particle volumes.
    equivalent_density = 0.010 / bed_volume
    if flow_basis == "main_text_120_Nml_min":
        molar_flow = 120.0e-6 / 60.0 * 101325.0 / (R * 273.15)
    elif flow_basis == "main_text_total_200_Nml_min":
        # The article text describes a two-component feed made by starting
        # N2 at 120 Nml/min and then adding CO2 at 80 Nml/min.  This explicit
        # basis preserves that reading as a diagnostic alternative; the
        # original 120 Nml/min basis remains the default historical adapter.
        molar_flow = 200.0e-6 / 60.0 * 101325.0 / (R * 273.15)
    elif flow_basis == "supplementary_484_mmol_h":
        molar_flow = 484.345157e-3 / 3600.0
    else:
        raise ValueError("Unknown Birmingham flow_basis")
    case = {
        "schema": "dynamic-psa-case/v1",
        "case_id": f"bham_2019_co2_n2_breakthrough_{n_cells}cell",
        "reporting_species": "CO2",
        "species": [
            {"name": "CO2", "mw_kg_mol": 0.0440095, "cp_j_mol_k": 37.13},
            {"name": "N2", "mw_kg_mol": 0.0280134, "cp_j_mol_k": 29.12},
        ],
        "bed": {
            "n_cells": int(n_cells), "length_m": 0.069, "diameter_m": 0.025,
            "void_fraction": 0.48, "pellet_void_fraction": 0.74,
            "bulk_density_kg_m3": equivalent_density, "particle_diameter_m": 0.001,
            "solid_cp_j_kg_k": 1000.0, "bed_k_w_m_k": 1.5,
            "axial_dispersion_m2_s": 5.0e-5,
            "wall_thickness_m": 0.002, "wall_density_kg_m3": 7700.0,
            "wall_cp_j_kg_k": 460.0, "wall_k_w_m_k": 10.9,
            "h_bed_wall_w_m2_k": 500000.0, "h_ambient_w_m2_k": 0.0,
            "ambient_k": tref, "kinetics": "loading_ldf",
            "energy_mode": "internal_energy", "materials": [material],
        },
        "feed": {"mole_fractions": [0.4, 0.6], "temperature_k": tref,
                 "pressure_bar": 25.0, "molar_flow_mol_s": molar_flow},
        "high_pressure_bar": 25.0, "low_pressure_bar": 1.0,
        "product_pressure_bar": 25.0, "initial_pressure_bar": 25.0,
        "initial_mole_fractions": [0.0, 1.0], "purge_to_feed_ratio": 0.0,
        "purge_gain_mol_s_bar": 0.0, "boundary_mode": "isobaric_breakthrough",
        "pressure_control_gain_mol_s_bar": 1.0e-3,
        "valve_regularization_bar": 1e-6,
        "valves": {"feed": {"taper_band_bar": 1.0}, "product": {"gain_mol_s_bar": 0.0},
                    "tail": {"gain_mol_s_bar": 0.0}, "equalization": {"gain_mol_s_sqrt_bar": 0.0}},
        "schedule": {"name": "breakthrough", "duration_s": 3600.0},
        "solver": {"rtol": 2e-7, "max_step_s": 8.0, "maximum_cycles": 1,
                    "minimum_cycles": 1, "css_tolerance": 1e-4,
                    "css_consecutive_cycles": 1, "mass_balance_tolerance": 1e-5,
                    "energy_balance_tolerance": 1e-5, "record_trace": True,
                    "trace_interval_s": 2.0},
        "provenance": {
            "source": "Birmingham UBIRe eData record 344",
            "doi": "10.25500/edata.bham.00000344",
            "article_doi": "10.1016/j.fuel.2019.05.095",
            "feed_basis": "40/60 mol% CO2/N2, 25 bar, 298.15 K",
            "flow_basis": flow_basis,
            "geometry_basis": "bed length 0.069 m and L/D 2.76; diameter derived as 0.025 m",
            "adsorbent_mass_basis": "measured 10 g charge converted to equivalent engine bulk-volume density; printed source particle-density field 262 kg/m3 retained as an audit value",
            "axial_dispersion_basis": "source glass-bead dispersion estimate approximately 5e-5 m2/s; source unit/exponent retained for audit and requires final table verification",
            "parameter_isolation": "CO2/N2 source branch; not reused for H2, CO2/CH4, CALF-20, or design cases",
            "n2_kinetic_status": "adapter_assumption_pending_source_resolution",
        },
        "claim": "Independent CO2/N2 isothermal breakthrough adapter; not yet a release-qualified source validation until N2 kinetic convention and exact thermal geometry are closed.",
    }
    parameter_payload = {"species": case["species"], "bed": case["bed"], "feed": case["feed"],
                         "initial_mole_fractions": case["initial_mole_fractions"],
                         "boundary_mode": case["boundary_mode"]}
    case["source_id"] = "birmingham_2019_co2_n2_ubire_344"
    case["adsorbent_id"] = "AC-MEA-MDEA"
    case["parameter_hash"] = hashlib.sha256(
        json.dumps(parameter_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return case


def bham_2015_mac_breakthrough_case(feed_co2=0.5, carrier="N2", n_cells=16,
                                    ldf_co2_per_s=1.0, axial_dispersion_m2_s=0.0):
    """Source-separated adapter for the Caldwell 2015 MAC breakthrough data.

    The released 2015 workbook contains MAC CO2/N2 and CO2/H2 traces at
    2,500 kPa, 298 K and 200 N mL min-1.  Tables S3/S4 provide the MAC
    dual-site equilibrium values at 25 degC, so those values are converted
    from Pa^-1 to bar^-1 and used directly at the source temperature.  H2 is
    intentionally represented as an inert carrier because the paper does not
    release matched pure-H2 equilibrium data.  The CO2 LDF coefficient is an
    explicit adapter assumption: the paper discusses small transport
    spreading but does not publish a fitted kinetic coefficient.

    This factory is diagnostic only.  It does not claim a source-equivalent
    kinetic replay, and its MAC parameters must not be copied into Wang,
    Otway, CALF-20 or NASA branches.
    """
    import math

    if carrier not in {"N2", "H2"}:
        raise ValueError("carrier must be N2 or H2")
    feed_co2 = float(feed_co2)
    if not 0.0 < feed_co2 < 1.0:
        raise ValueError("feed_co2 must lie strictly between 0 and 1")
    if int(n_cells) < 2:
        raise ValueError("n_cells must be at least 2")
    if ldf_co2_per_s <= 0.0:
        raise ValueError("ldf_co2_per_s must be positive")
    # Table S4 at 25 degC, converted from Pa^-1 to bar^-1.  Binding energies
    # are set to zero because this branch is deliberately isothermal; the
    # source heat-of-adsorption values are retained in provenance only.
    co2_a1 = [0.81, 8.77]
    co2_b0_bar = [1.17e-4 * 1.0e5, 3.00e-6 * 1.0e5]
    if carrier == "N2":
        carrier_a1 = [3.74, 172.0]
        carrier_b0_bar = [1.06e-6 * 1.0e5, 1.04e-9 * 1.0e5]
    else:
        # The source has CO2/H2 breakthrough curves but no matched H2
        # equilibrium table.  Zero capacity is an explicit inert-carrier
        # closure, not an invented H2 adsorption fit.
        carrier_a1 = [0.0, 0.0]
        carrier_b0_bar = [0.0, 0.0]
    names = ["CO2", carrier]
    carrier_mw = 0.0280134 if carrier == "N2" else 0.00201588
    carrier_cp = 29.12 if carrier == "N2" else 28.84
    material = {
        "name": "Caldwell2015_MAC_25C_source_equivalent_equilibrium",
        "end_fraction": 1.0,
        "a1": [co2_a1, carrier_a1],
        "a2": [[0.0, 0.0], [0.0, 0.0]],
        "b0_per_bar": [co2_b0_bar, carrier_b0_bar],
        "binding_energy_j_mol": [[0.0, 0.0], [0.0, 0.0]],
        "ldf_per_s": [float(ldf_co2_per_s), 0.0],
        "source_parameter_status": "equilibrium_source_equivalent;_CO2_LDF_and_H2_or_N2_transport_are_adapter_closures",
    }
    molar_flow = 200.0e-6 / 60.0 * 101325.0 / (R * 273.15)
    case = {
        "schema": "dynamic-psa-case/v1",
        "case_id": f"bham_2015_mac_co2_{carrier.lower()}_breakthrough_{feed_co2:.3f}_{int(n_cells)}cell",
        "source_id": "birmingham_2015_co2_n2_h2_figure_data",
        "adsorbent_id": "Modified Activated Carbon",
        "reporting_species": "CO2",
        "species": [
            {"name": "CO2", "mw_kg_mol": 0.0440095, "cp_j_mol_k": 37.13},
            {"name": carrier, "mw_kg_mol": carrier_mw, "cp_j_mol_k": carrier_cp},
        ],
        "bed": {
            "n_cells": int(n_cells), "length_m": 0.065, "diameter_m": 0.025,
            # The article reports bed voidage but not intrapellet voidage.
            # The positive value below is an explicit transport-adapter
            # closure required by the engine's partial-pressure LDF state.
            "void_fraction": 0.40, "pellet_void_fraction": 0.74,
            "bulk_density_kg_m3": 191.0, "particle_diameter_m": 0.001,
            "solid_cp_j_kg_k": 1000.0, "bed_k_w_m_k": 1.5,
            "axial_dispersion_m2_s": float(axial_dispersion_m2_s),
            "wall_thickness_m": 0.002, "wall_density_kg_m3": 7700.0,
            "wall_cp_j_kg_k": 460.0, "wall_k_w_m_k": 10.9,
            "h_bed_wall_w_m2_k": 500000.0, "h_ambient_w_m2_k": 500000.0,
            "ambient_k": 298.15, "kinetics": "partial_pressure_ldf",
            "energy_mode": "source_sensible", "materials": [material],
        },
        "feed": {"mole_fractions": [feed_co2, 1.0 - feed_co2],
                 "temperature_k": 298.15, "pressure_bar": 25.0,
                 "molar_flow_mol_s": molar_flow},
        "high_pressure_bar": 25.0, "low_pressure_bar": 1.0,
        "product_pressure_bar": 25.0, "initial_pressure_bar": 25.0,
        "initial_mole_fractions": [0.0, 1.0], "purge_to_feed_ratio": 0.0,
        "purge_gain_mol_s_bar": 0.0, "boundary_mode": "isobaric_breakthrough",
        "pressure_control_gain_mol_s_bar": 1.0e-3,
        "valve_regularization_bar": 1e-6,
        "valves": {"feed": {"taper_band_bar": 1.0}, "product": {"gain_mol_s_bar": 0.0},
                    "tail": {"gain_mol_s_bar": 0.0}, "equalization": {"gain_mol_s_sqrt_bar": 0.0}},
        "schedule": {"name": "breakthrough", "duration_s": 5400.0},
        "solver": {"rtol": 2e-7, "max_step_s": 8.0, "maximum_cycles": 1,
                    "minimum_cycles": 1, "css_tolerance": 1e-4,
                    "css_consecutive_cycles": 1, "mass_balance_tolerance": 1e-5,
                    "energy_balance_tolerance": 1e-5, "record_trace": True,
                    "trace_interval_s": 10.0},
        "provenance": {
            "article_doi": "10.1021/acs.energyfuels.5b00164",
            "data_doi": "10.25500/eData.bham.00000060",
            "equilibrium_basis": "Tables S3/S4, MAC 25 degC; B in Pa^-1 converted to bar^-1",
            "geometry_basis": "Table 5: 0.065 m bed, 0.025 m diameter, 6.11 g MAC, voidage 0.40",
            "feed_basis": "Table 5: 2500 kPa, 298 K, 200 N mL min^-1",
            "heat_of_adsorption_source_kj_mol": {"CO2": 26.6, "N2": 15.7},
            "isothermal_closure": "source_sensible at 298.15 K with zero binding-energy source term",
            "kinetic_status": "unknown_source_LDF; explicit adapter value",
            "pellet_voidage_status": "not_reported_in_source; 0.74 is an adapter closure required for the partial-pressure LDF state",
            "carrier_status": "N2 uses Table S3 equilibrium; H2 is inert because matched H2 equilibrium is not reported",
            "parameter_isolation": "Birmingham 2015 CCS branch only; never reused by Wang, Otway, CALF-20 or NASA",
        },
        "claim": "Diagnostic source-specific MAC equilibrium adapter; not a release-qualified transport or non-isothermal validation.",
    }
    payload = {"species": case["species"], "bed": case["bed"], "feed": case["feed"],
               "boundary_mode": case["boundary_mode"], "source_id": case["source_id"]}
    case["parameter_hash"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return case


def tamof1_fixed_bed_case(*, pressure_bar: float, temperature_k: float,
                          molar_flow_mol_s: float, ldf_per_s,
                          particle_diameter_m: float, duration_s: float,
                          bed_mass_g: float = 0.7, bed_length_m: float = 0.01,
                          n_cells: int = 16, void_fraction: float = 0.30,
                          pellet_void_fraction: float = 0.598,
                          feed_mole_fractions=(0.50, 0.50, 0.0),
                          initial_mole_fractions=(0.0, 0.0, 1.0),
                          kinetics: str = "partial_pressure_ldf"):
    """Build an explicitly parameterized TAMOF-1 single-bed diagnostic case.

    The Nature Communications article publishes the two bed masses/lengths,
    dehydrated bulk density and material porosity, but not a source-fitted LDF
    coefficient or a complete dynamic pressure/thermal closure.  Therefore
    every dynamic adapter value that is not directly published is required as
    a keyword argument; no hidden default is supplied.  The factory is
    isothermal and fail-closed outside the 303 K equilibrium namespace.  Its
    result is diagnostic-only and cannot enter cyclic CSS, DWSIM, TEA/LCA or
    training release.
    """
    import math
    from .tamof1_2025_adapter import tamof1_303k_material

    pressure_bar = float(pressure_bar)
    temperature_k = float(temperature_k)
    molar_flow_mol_s = float(molar_flow_mol_s)
    particle_diameter_m = float(particle_diameter_m)
    duration_s = float(duration_s)
    bed_mass_g = float(bed_mass_g)
    bed_length_m = float(bed_length_m)
    void_fraction = float(void_fraction)
    pellet_void_fraction = float(pellet_void_fraction)
    feed_mole_fractions = [float(value) for value in feed_mole_fractions]
    initial_mole_fractions = [float(value) for value in initial_mole_fractions]
    if abs(temperature_k - 303.0) > 1.0e-9:
        raise ValueError("TAMOF-1 dynamic case is closed outside the published 303 K equilibrium namespace")
    if kinetics not in {"partial_pressure_ldf", "loading_ldf"}:
        raise ValueError("TAMOF-1 diagnostic kinetics must be partial_pressure_ldf or loading_ldf")
    if pressure_bar <= 0 or molar_flow_mol_s <= 0 or particle_diameter_m <= 0 or duration_s <= 0:
        raise ValueError("pressure, flow, particle diameter and duration must be positive")
    if bed_mass_g <= 0 or bed_length_m <= 0 or not 0 < void_fraction < 1 or not 0 < pellet_void_fraction < 1:
        raise ValueError("invalid TAMOF-1 bed geometry or porosity")
    ldf = [float(value) for value in ldf_per_s]
    if len(ldf) != 3 or any(value <= 0 for value in ldf):
        raise ValueError("ldf_per_s must contain three positive, explicitly supplied values for CO2/CH4/N2")
    for label, composition in (("feed_mole_fractions", feed_mole_fractions),
                               ("initial_mole_fractions", initial_mole_fractions)):
        if len(composition) != 3 or any(value < 0.0 for value in composition) or not math.isclose(sum(composition), 1.0, rel_tol=0.0, abs_tol=1.0e-10):
            raise ValueError(f"{label} must contain three non-negative mole fractions summing to one")
    # The article reports dehydrated bulk density = 0.530 g cm^-3.  The
    # engine's bulk-density field is multiplied by (1-bed_void_fraction), so
    # derive the cross-section from the published mass and length rather than
    # inserting an unrelated tube diameter.
    source_bulk_density_kg_m3 = 530.0
    # Bed.mass in the common engine is bulk_density * total cell volume,
    # while the article's 0.530 g cm^-3 is the dehydrated packed-material
    # density before applying the explicitly declared bed void fraction.
    # Convert once here so the reconstructed bed mass remains 0.7/7.7 g.
    engine_bulk_density_kg_m3 = source_bulk_density_kg_m3 * (1.0 - void_fraction)
    mass_kg = bed_mass_g * 1.0e-3
    area = mass_kg / (engine_bulk_density_kg_m3 * bed_length_m)
    diameter_m = math.sqrt(4.0 * area / math.pi)
    material = tamof1_303k_material()
    material.update({
        "name": "TAMOF1_2025_Table16_dynamic_diagnostic",
        "end_fraction": 1.0,
        "ldf_per_s": ldf,
        # No source-equivalent heat-of-adsorption/thermal storage mapping is
        # published for this dynamic replay. Zero values keep the isothermal
        # diagnostic explicit and prevent accidental non-isothermal claims.
        "binding_energy_j_mol": ([[0.0] for _ in range(3)] if kinetics == "loading_ldf"
                                   else [[0.0, 0.0] for _ in range(3)]),
        "source_parameter_status": "303_K_equilibrium_source;_LDF_and_dynamic_closures_explicit_adapter_inputs",
    })
    # Case JSON is the restart/provenance boundary; do not leave NumPy arrays
    # inside a serialized input even though the equilibrium adapter accepts
    # them in memory.
    for key in ("tamof_qm_mol_kg", "tamof_k_per_bar_power", "tamof_n"):
        material[key] = material[key].tolist()
    species = [
        {"name": "CO2", "mw_kg_mol": 0.0440095, "cp_j_mol_k": 37.13},
        {"name": "CH4", "mw_kg_mol": 0.0160425, "cp_j_mol_k": 35.69},
        {"name": "N2", "mw_kg_mol": 0.0280134, "cp_j_mol_k": 29.12},
    ]
    case = {
        "schema": "dynamic-psa-case/v1",
        "case_id": f"tamof1_303k_breakthrough_{bed_mass_g:g}g_{int(n_cells)}cell",
        "source_id": "tamof1_2025_co2_ch4_family_diagnostic",
        "adsorbent_id": "TAMOF-1",
        "reporting_species": "CO2",
        "species": species,
        "bed": {
            "n_cells": int(n_cells), "length_m": bed_length_m, "diameter_m": diameter_m,
            "void_fraction": void_fraction, "pellet_void_fraction": pellet_void_fraction,
            "bulk_density_kg_m3": engine_bulk_density_kg_m3, "particle_diameter_m": particle_diameter_m,
            "solid_cp_j_kg_k": 1000.0, "bed_k_w_m_k": 0.2,
            "wall_thickness_m": 0.001, "wall_density_kg_m3": 7850.0,
            "wall_cp_j_kg_k": 460.0, "wall_k_w_m_k": 10.9,
            "h_bed_wall_w_m2_k": 0.0, "h_ambient_w_m2_k": 0.0,
            "ambient_k": 303.0,
            "energy_mode": "source_sensible", "materials": [material],
            "kinetics": kinetics,
        },
        "feed": {"mole_fractions": feed_mole_fractions, "temperature_k": temperature_k,
                 "pressure_bar": pressure_bar, "molar_flow_mol_s": molar_flow_mol_s},
        "high_pressure_bar": pressure_bar, "low_pressure_bar": 0.01,
        "product_pressure_bar": pressure_bar, "initial_pressure_bar": pressure_bar,
        "initial_mole_fractions": initial_mole_fractions,
        "purge_to_feed_ratio": 0.0, "purge_gain_mol_s_bar": 0.0,
        "boundary_mode": "isobaric_breakthrough",
        "pressure_control_gain_mol_s_bar": 1.0e-3,
        "valve_regularization_bar": 1.0e-6,
        "valves": {"feed": {"taper_band_bar": 0.05}, "product": {"gain_mol_s_bar": 0.0},
                    "tail": {"gain_mol_s_bar": 0.0}, "equalization": {"gain_mol_s_sqrt_bar": 0.0}},
        "schedule": {"name": "breakthrough", "duration_s": duration_s},
        "solver": {"rtol": 2e-7, "max_step_s": 8.0, "maximum_cycles": 1,
                    "minimum_cycles": 1, "css_tolerance": 1e-4,
                    "css_consecutive_cycles": 1, "mass_balance_tolerance": 1e-5,
                    "energy_balance_tolerance": 1e-5, "record_trace": True,
                    "trace_interval_s": 5.0},
        "provenance": {
            "article_doi": "10.1038/s41467-025-58426-w",
            "data_doi": "10.5281/zenodo.13126536",
            "geometry_basis": "Published column mass/bed length and dehydrated bulk density; engine bulk-density basis converted once with explicit bed void fraction",
            "source_dehydrated_bulk_density_kg_m3": source_bulk_density_kg_m3,
            "bed_mass_g": bed_mass_g,
            "material_porosity_basis": "Published TAMOF-1 porosity 0.598 used as pellet_void_fraction",
            "equilibrium_basis": "Supplementary Table 16, 303 K dual-site Langmuir-Freundlich",
            "dynamic_inputs": {"ldf_per_s": ldf, "particle_diameter_m": particle_diameter_m, "molar_flow_mol_s": molar_flow_mol_s, "feed_mole_fractions": feed_mole_fractions, "initial_mole_fractions": initial_mole_fractions},
            "boundary_basis": "Isobaric breakthrough diagnostic; source pressure/flow boundary mapping remains open",
            "thermal_basis": "Isothermal source_sensible closure with zero binding-energy term; not a non-isothermal validation",
            "parameter_isolation": "TAMOF1_CO2_CH4 only; never mix with Wang, CCS, CALF-20 or other biogas branches",
        },
        "claim": "TAMOF-1 303 K single-bed diagnostic case; explicit LDF/geometry closures; not source-equivalent validation or cyclic design release.",
        "release_qualified": False,
    }
    payload = {"species": case["species"], "bed": case["bed"], "feed": case["feed"],
               "initial_mole_fractions": case["initial_mole_fractions"],
               "boundary_mode": case["boundary_mode"], "source_id": case["source_id"]}
    case["parameter_hash"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"),
                   default=lambda value: value.tolist() if hasattr(value, "tolist") else str(value)).encode("utf-8")
    ).hexdigest()
    return case
