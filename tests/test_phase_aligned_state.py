import pytest

from psa_dynamic.phase_aligned_state import resolve_phase_aligned_states


def _samples(pressure, temperature):
    return [
        {"time_s": 0.0, "duration_s": 1.0, "flow_mol_s": 2.0,
         "pressure_bar": pressure, "temperature_k": temperature},
        {"time_s": 1.0, "duration_s": 1.0, "flow_mol_s": 1.0,
         "pressure_bar": pressure + 1.0, "temperature_k": temperature + 3.0},
    ]


def test_phase_aligned_states_are_molar_exposure_weighted():
    result = {"phase_aligned_streams": {
        "FEED": _samples(10.0, 300.0),
        "PRODUCT": _samples(9.0, 301.0),
        "TAIL": _samples(1.0, 299.0),
    }}
    report = resolve_phase_aligned_states(result)
    assert report["status"] == "resolved"
    assert report["streams"]["PRODUCT"]["pressure_bar"] == pytest.approx(9.0 + 1.0 / 3.0)
    assert report["streams"]["PRODUCT"]["temperature_k"] == pytest.approx(302.0)


def test_missing_or_invalid_streams_remain_unresolved():
    result = {"phase_aligned_streams": {
        "FEED": _samples(10.0, 300.0),
        "PRODUCT": [{"time_s": 1.0, "duration_s": 1.0, "flow_mol_s": 1.0,
                     "pressure_bar": 9.0, "temperature_k": 301.0},
                    {"time_s": 0.0, "duration_s": 1.0, "flow_mol_s": 1.0,
                     "pressure_bar": 9.0, "temperature_k": 301.0}],
    }}
    report = resolve_phase_aligned_states(result)
    assert report["status"] == "unresolved"
    assert "samples_missing_or_empty" in report["streams"]["TAIL"]["errors"]
    assert "sample_1:nonmonotonic_time_s" in report["streams"]["PRODUCT"]["errors"]
