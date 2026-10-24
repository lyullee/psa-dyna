import pytest

from psa_dynamic.schedules import build_schedule, validate_schedule_contract


def test_guderian_six_named_operations_are_two_beds_on_four_global_intervals():
    schedule = build_schedule({
        "name": "guderian_two_bed_6step_source",
        "half_cycle_s": 60.0,
        "equalization_s": 1.0,
    })
    contract = validate_schedule_contract(schedule)
    assert contract["counts_are_independent"] is True
    assert schedule["number_of_beds"] == 2
    assert schedule["number_of_cycle_steps"] == 4
    assert schedule["source_named_operations_per_bed"] == 6
    contract = validate_schedule_contract(schedule)
    assert contract["number_of_beds"] == 2
    assert contract["number_of_cycle_steps"] == 4
    assert contract["source_named_operations_per_bed"] == 6
    assert contract["source_operation_count_definition"] == "distinct_named_operations_executed_by_each_bed"
    assert schedule["cycle_time_s"] == 120.0
    assert [slot["duration_s"] for slot in schedule["slots"]] == [1.0, 59.0, 1.0, 59.0]
    assert [slot["actions"] for slot in schedule["slots"]] == [
        ["EQ-D", "EQ-R"], ["BR", "AP"],
        ["EQ-R", "EQ-D"], ["AP", "BR"],
    ]
    # The high-pressure equalization donor must regenerate next, while the
    # low-pressure receiver must adsorb/produce.  This is the phase relation
    # visible in Marcinek et al. Figure 7 and Figure 8.
    for equalization_index, process_index in ((0, 1), (2, 3)):
        connections = schedule["slots"][equalization_index]["equalization_connections"]
        donor = connections[0]["donor"]
        receiver = connections[0]["receiver"]
        actions = schedule["slots"][process_index]["actions"]
        assert actions[donor] == "BR"
        assert actions[receiver] == "AP"


def test_guderian_equalization_retains_both_source_connections():
    schedule = build_schedule({
        "name": "guderian_two_bed_6step_source",
        "half_cycle_s": 45.0,
        "equalization_s": 1.0,
    })
    first = schedule["slots"][0]["equalization_connections"]
    second = schedule["slots"][2]["equalization_connections"]
    assert {(row["donor_endpoint"], row["receiver_endpoint"]) for row in first} == {
        ("top", "top"), ("bottom", "bottom"),
    }
    assert {(row["donor"], row["receiver"]) for row in first} == {(0, 1)}
    assert {(row["donor"], row["receiver"]) for row in second} == {(1, 0)}
    assert schedule["source_equivalent"] is False


def test_guderian_schedule_fails_closed_on_impossible_timing():
    with pytest.raises(ValueError, match="must exceed"):
        build_schedule({
            "name": "guderian_two_bed_6step_source",
            "half_cycle_s": 1.0,
            "equalization_s": 1.0,
        })
