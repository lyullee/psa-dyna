"""Source-local contract for Dantas et al. (2011) CO2/N2--zeolite 13X.

The SciELO article reports fixed-bed breakthrough experiments at several
temperatures, the Toth equilibrium parameters, LDF coefficients, bed
properties, and the energy/momentum model assumptions.  It is therefore a
useful independent CCS *single-bed* branch.  It is not a cyclic VPSA
historian, and no value in this module may be copied into Wang H2,
Edinburgh/biogas, CALF-20, or another CCS namespace.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


SOURCE_ID = "DANTAS_2011_CO2_N2_ZEOLITE_13X"
DOI = "10.1590/S0104-66322011000300018"
PUBLIC_URL = "https://www.scielo.br/j/bjce/a/RBxBDJ3CVLbYJSmrzwPw9sN/"
# The source contains both CO2/He carrier-gas runs (1--4) and CO2/N2 runs
# (5--12).  Helium is inert in the source model and has no Toth parameter;
# it is still retained in the component contract so the carrier gas cannot be
# silently confused with nitrogen.
COMPONENTS = ("CO2", "N2", "He")
ARTICLE_PATH = Path(__file__).resolve().parent / "literature" / "source_dantas_2011_co2_n2_13x.pdf"
ARTICLE_SHA256 = "F1D580829CA31D3F5F224BAA8FE8B645BE27C7686C04922F0B0F1780F20C3F4E"


class DantasSourceError(ValueError):
    """Raised when the Dantas source is used outside its published contract."""


# Table 6.  K0 is the source's bar^-1 prefactor and heat is stored as the
# positive magnitude of -Delta H, in kJ/mol.  The source applies the van't
# Hoff temperature dependence to K; no parameter is inferred from another
# publication.
_TOTH = {
    "CO2": {"qm_mol_kg": 5.09, "n": 0.429, "K0_bar_inv": 4.31e-4, "heat_kJ_mol": 29.38},
    "N2": {"qm_mol_kg": 3.08, "n": 0.869, "K0_bar_inv": 8.81e-5, "heat_kJ_mol": 17.19},
}

_LDF_S_INV = {
    1: {"CO2": 0.226},
    2: {"CO2": 0.370},
    3: {"CO2": 1.196},
    4: {"CO2": 2.525},
    5: {"CO2": 0.065, "N2": 2.733},
    6: {"CO2": 0.108, "N2": 5.102},
    7: {"CO2": 0.353, "N2": 9.566},
    8: {"CO2": 0.758, "N2": 31.844},
    9: {"CO2": 0.067, "N2": 2.792},
    10: {"CO2": 0.106, "N2": 5.014},
    11: {"CO2": 0.344, "N2": 9.324},
    12: {"CO2": 0.726, "N2": 30.519},
}

_CASE_TEMPERATURE_K = {1: 301.0, 2: 323.0, 3: 373.0, 4: 423.0,
                       5: 301.0, 6: 323.0, 7: 373.0, 8: 423.0,
                       9: 306.0, 10: 323.0, 11: 373.0, 12: 423.0}


@dataclass(frozen=True)
class DantasCase:
    """One published Dantas fixed-bed run, without cross-source defaults."""

    run: int

    def manifest(self) -> dict[str, Any]:
        if self.run not in _CASE_TEMPERATURE_K:
            raise DantasSourceError("Dantas run must be one of 1--12")
        binary = self.run >= 5
        nitrogen_saturated = self.run >= 9
        return {
            "run": self.run,
            "temperature_K": _CASE_TEMPERATURE_K[self.run],
            "pressure_bar": 1.20 if nitrogen_saturated else 1.02,
            "feed_mole_fraction": {"CO2": 0.10 if nitrogen_saturated else 0.20,
                                    "N2": 0.90 if binary else 0.0,
                                    "He": 0.80 if not binary else 0.0},
            "flow_m3_s": 3.5e-5 if nitrogen_saturated else (1.0e-6 if binary else 5.0e-7),
            "energy_boundary": "non-adiabatic" if nitrogen_saturated else "adiabatic",
            "initial_condition": "N2-saturated bed" if nitrogen_saturated else "He-regenerated bed",
            "ldf_overall_s_inv": dict(_LDF_S_INV[self.run]),
            "source_parameter_namespace": SOURCE_ID,
        }


@dataclass(frozen=True)
class DantasCO2N2Zeolite13X:
    """Hash-pinned Dantas source bundle for independent CCS diagnostics."""

    gas_components: tuple[str, ...] = COMPONENTS

    def _check(self) -> None:
        if self.gas_components != COMPONENTS:
            raise DantasSourceError("Dantas source is restricted to CO2/N2/He component order")
        if not ARTICLE_PATH.is_file():
            raise DantasSourceError(f"Dantas source file missing: {ARTICLE_PATH}")
        digest = hashlib.sha256(ARTICLE_PATH.read_bytes()).hexdigest().upper()
        if digest != ARTICLE_SHA256:
            raise DantasSourceError(f"Dantas source SHA-256 mismatch: {ARTICLE_PATH.name}")

    def toth_parameters(self) -> dict[str, dict[str, float]]:
        self._check()
        return {name: dict(values) for name, values in _TOTH.items()}

    def toth_loading_mol_kg(self, component: str, pressure_bar: float, temperature_K: float) -> float:
        """Evaluate the source Toth form using only Dantas parameters."""
        params = self.toth_parameters().get(component)
        if params is None:
            raise DantasSourceError(f"unknown Dantas component: {component}")
        if pressure_bar < 0.0 or temperature_K <= 0.0:
            raise DantasSourceError("pressure must be non-negative and temperature positive")
        R = 8.31446261815324  # J mol^-1 K^-1
        # The paper defines K=K0*exp(-DeltaH/(R*T)); Table 6 reports the
        # positive magnitude ``-DeltaH``. Therefore use exp(+|DeltaH|/(R*T)).
        K = params["K0_bar_inv"] * math.exp(params["heat_kJ_mol"] * 1000.0 / (R * temperature_K))
        kp = K * pressure_bar
        n = params["n"]
        return params["qm_mol_kg"] * kp / (1.0 + kp**n) ** (1.0 / n)

    def bed_parameters(self, nitrogen_saturated: bool = False) -> dict[str, float]:
        self._check()
        if nitrogen_saturated:
            return {
                "bed_length_m": 0.83, "bed_diameter_m": 0.021,
                "bed_mass_kg": 0.210, "bed_porosity": 0.41,
                "wall_thickness_m": 0.0041, "wall_heat_capacity_J_kgK": 500.0,
                "wall_conductivity_W_mK": 13.4, "wall_density_kg_m3": 8238.0,
            }
        return {
            "bed_length_m": 0.171, "bed_diameter_m": 0.022,
            "bed_mass_kg": 0.0455, "solid_density_kg_m3": 1950.0,
            "particle_density_kg_m3": 1228.5, "particle_diameter_m": 0.0029,
            "particle_porosity": 0.37, "bed_porosity": 0.43,
            "tortuosity": 2.2, "solid_heat_capacity_J_kgK": 920.0,
            "wall_thickness_m": 0.0015, "wall_heat_capacity_J_kgK": 440.0,
            "wall_density_kg_m3": 7830.0,
        }

    def manifest(self) -> dict[str, Any]:
        self._check()
        return {
            "schema": "dynamic-psa-dantas-2011-co2-n2-source-adapter/v1",
            "source_id": SOURCE_ID, "doi": DOI, "public_url": PUBLIC_URL,
            "article_sha256": ARTICLE_SHA256,
            "gas_components": list(self.gas_components),
            "equilibrium": {"model": "Toth + van't Hoff", "parameters": self.toth_parameters()},
            "runs": [DantasCase(run).manifest() for run in range(1, 13)],
            "model": {"mass_transfer": "LDF", "energy": "gas/solid/wall; adiabatic or non-adiabatic by run", "momentum": "included for mixtures", "gas_phase": "ideal gas"},
            "release_flags": {
                "source_equilibrium_released": True,
                "source_single_bed_breakthrough_released": False,
                "source_cyclic_vpsa_released": False,
                "raw_machine_readable_historian": False,
                "dwsim_application_allowed": False,
                "tea_lca_allowed": False,
                "training_eligible": False,
            },
            "parameter_isolation": "DANTAS_2011_CO2_N2_13X_ONLY_no_Wang_Edinburgh_Buchner_CALF20_NASA_transfer",
            "claim_boundary": "Dantas provides source-specific fixed-bed CO2/N2--13X conditions and model parameters. Experimental curves are published as figures; no raw historian or cyclic CSS validation is released.",
        }

    def build_engine_material(self, components: tuple[str, ...] = ("CO2", "N2", "He")) -> dict[str, Any]:
        """Return the Dantas Toth material for the common engine.

        This is an engine-material contract, not an experimental-release
        claim.  The source uses pure-component Toth fits; the carrier gas in
        each run is represented by an explicit zero-capacity inert row.
        """
        self._check()
        requested = tuple(components)
        if requested not in (("CO2", "N2", "He"), ("CO2", "N2"), ("CO2", "He")):
            raise DantasSourceError("Dantas engine material supports only CO2/N2/He, CO2/N2 or CO2/He order")
        values = {name: _TOTH.get(name, {"qm_mol_kg": 0.0, "n": 1.0, "K0_bar_inv": 0.0, "heat_kJ_mol": 0.0}) for name in requested}
        return {
            "name": "Dantas2011_Zeolite13X_Toth_source_local",
            "model": "dantas_toth_source",
            "end_fraction": 1.0,
            "dantas_qm_mol_kg": [[float(values[name]["qm_mol_kg"])] for name in requested],
            "dantas_n": [[float(values[name]["n"])] for name in requested],
            "dantas_k0_per_bar": [[float(values[name]["K0_bar_inv"])] for name in requested],
            "dantas_heat_j_mol": [[float(values[name]["heat_kJ_mol"] * 1000.0)] for name in requested],
            # The loading-LDF energy ledger uses the common binding-energy
            # field.  It is copied only within this source-local material;
            # no Dantas heat is reused by another gas/adsorbent branch.
            "binding_energy_j_mol": [[float(values[name]["heat_kJ_mol"] * 1000.0)] for name in requested],
            "parameter_provenance": SOURCE_ID,
        }


def source_manifest() -> dict[str, Any]:
    return DantasCO2N2Zeolite13X().manifest()
