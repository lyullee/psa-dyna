import pytest

from psa_dynamic.cases import wang_case
from psa_dynamic.core import Network
from psa_dynamic.runner import atomic_json, run_case
from psa_dynamic.schedules import build_schedule, validate_schedule_contract


def test_bed_count_and_cycle_step_count_are_independent_for_4_and_8_step_plugins():
    four = build_schedule({
        "name": "skarstrom_4step",
        "pressurization_s": 10.0,
        "adsorption_s": 30.0,
    })
    eight = build_schedule({"name": "equalized_8step", "slot_s": 10.0})

    assert (four["number_of_beds"], four["number_of_cycle_steps"]) == (2, 4)
    assert (eight["number_of_beds"], eight["number_of_cycle_steps"]) == (4, 8)
    for schedule in (four, eight):
        contract = validate_schedule_contract(schedule)
        assert contract == {
            "number_of_beds": schedule["number_of_beds"],
            "number_of_cycle_steps": schedule["number_of_cycle_steps"],
            "bed_count_definition": "action_rows",
            "cycle_step_count_definition": "slot_count",
            "counts_are_independent": True,
        }
        assert len(schedule["slots"]) == schedule["number_of_cycle_steps"]
        assert all(len(slot["actions"]) == schedule["number_of_beds"]
                   for slot in schedule["slots"])


def test_network_rejects_cycle_step_metadata_that_does_not_match_slots():
    case = wang_case("closed", 2)
    schedule = build_schedule(case["schedule"])
    schedule["number_of_cycle_steps"] = 2
    with pytest.raises(ValueError, match="cycle-step count mismatch"):
        Network(case, schedule)


def test_network_rejects_a_slot_with_the_wrong_number_of_bed_actions():
    case = wang_case("closed", 2)
    schedule = build_schedule(case["schedule"])
    schedule["number_of_beds"] = 2
    schedule["slots"][0]["actions"] = ["IDLE"]
    with pytest.raises(ValueError, match="bed count mismatch"):
        Network(case, schedule)


def test_network_rejects_cycle_time_that_is_not_slot_sum():
    case = wang_case("skarstrom_4step", 2)
    case["schedule"].update({"pressurization_s": 1.0, "adsorption_s": 2.0})
    schedule = build_schedule(case["schedule"])
    schedule["cycle_time_s"] += 1.0
    with pytest.raises(ValueError, match="cycle_time_s"):
        Network(case, schedule)


def test_executed_result_records_the_same_independent_count_definitions(tmp_path):
    case = wang_case("closed", 2)
    case_path = tmp_path / "case.json"
    atomic_json(case_path, case)
    result = run_case(case_path, tmp_path / "result")
    assert result["number_of_beds"] == 1
    assert result["number_of_cycle_steps"] == 1
    assert result["schedule_contract"] == {
        "number_of_beds": 1,
        "number_of_cycle_steps": 1,
        "bed_count_definition": "action_rows",
        "cycle_step_count_definition": "slot_count",
        "counts_are_independent": True,
    }
