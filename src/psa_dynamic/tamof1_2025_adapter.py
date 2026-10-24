"""Source-local TAMOF-1 equilibrium parameters.

The public TAMOF-1 supplementary information reports a dual-site
Langmuir--Freundlich fit at 303 K.  This adapter exposes only that published
equilibrium namespace.  It deliberately does not invent an LDF coefficient,
wall heat-transfer coefficient, pressure-drop law, or temperature dependence
for CH4/N2, so it cannot silently enter a non-isothermal PSA design case.
"""

from __future__ import annotations

import numpy as np


TAMOF1_ARTICLE_DOI = "10.1038/s41467-025-58426-w"
TAMOF1_DATA_DOI = "10.5281/zenodo.13126536"
TAMOF1_PARAMETER_TABLE = "Supplementary Table 16"
TAMOF1_REFERENCE_TEMPERATURE_K = 303.0

# Species order is deliberately explicit and source-local.
TAMOF1_SPECIES = ("CO2", "CH4", "N2")

# Table 16 values: qm [mmol g^-1 == mol kg^-1], K [bar^(-n)], n [-].
_QM = np.array([
    [1.34281, 4.64294],
    [0.41023, 0.86851],
    [0.17408, 0.73201],
], dtype=float)
_K = np.array([
    [56.61094, 1.02766],
    [3.73563, 0.74271],
    [0.27453, 0.39096],
], dtype=float)
_N = np.array([
    [1.12067, 0.82597],
    [0.82172, 0.46732],
    [0.40343, 0.96000],
], dtype=float)


def tamof1_303k_material() -> dict:
    """Return a fresh, immutable-by-convention TAMOF-1 equilibrium mapping."""
    return {
        "model": "tamof1_dual_site_langmuir_freundlich_303k",
        "species_order": list(TAMOF1_SPECIES),
        "tamof_qm_mol_kg": _QM.copy(),
        "tamof_k_per_bar_power": _K.copy(),
        "tamof_n": _N.copy(),
        "tamof_reference_temperature_k": TAMOF1_REFERENCE_TEMPERATURE_K,
        "tamof_temperature_tolerance_k": 1.0e-9,
        "parameter_provenance": (
            f"TAMOF-1; {TAMOF1_ARTICLE_DOI}; {TAMOF1_PARAMETER_TABLE}; "
            f"data DOI {TAMOF1_DATA_DOI}; 303 K"
        ),
        "release_qualified": False,
        "dynamic_use_allowed": False,
    }


__all__ = [
    "TAMOF1_ARTICLE_DOI",
    "TAMOF1_DATA_DOI",
    "TAMOF1_PARAMETER_TABLE",
    "TAMOF1_REFERENCE_TEMPERATURE_K",
    "TAMOF1_SPECIES",
    "tamof1_303k_material",
]
