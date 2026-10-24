"""Source-labelled two-column VPSA timing reconstruction.

Figure 3 of Streb & Mazzotti (2021, part 2) draws one half-cycle.  Eight
named operations occur *per bed across the full cycle*: PE-Press, Feed, Idle,
PE-BD, BD, HP, BD-vac, LP.  Feed includes pressurization and adsorption; its
internal split is not reported in Fig. 3.  Global event intervals are the
union of both columns' transition times and need not number eight.

This object is not itself a source-equivalent dynamic input.  A separate
diagnostic schedule adapter can execute its topology with explicit laboratory
bottle inlets, but the paper's pressure trajectory, Wilson/RAST equilibrium,
thermal closure and experimental comparison remain unqualified.  Table 5
rows 1-8 were visually checked against the original full-size table.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


SOURCE = Path(__file__).resolve().parent / "literature" / "STREB_MAZZOTTI_2021_SOURCE_CONTRACT_SCREEN.json"
NAMED_SEQUENCE = ("PE-Press", "Feed", "Idle", "PE-BD", "BD", "HP", "BD-vac", "LP")


def _operation_at(t: float, boundaries: tuple[float, ...], names: tuple[str, ...]) -> str:
    for i, name in enumerate(names):
        if boundaries[i] <= t < boundaries[i + 1]:
            return name
    raise ArithmeticError("Time lies outside the source half-cycle")


def _collapsed_actions(slots: list[dict], bed: int) -> list[str]:
    result = []
    for slot in slots:
        action = slot["actions"][bed]
        if not result or result[-1] != action:
            result.append(action)
    return result


def build_streb_two_bed_timeline(experiment: int = 1) -> dict:
    """Return source-checked timing for Table-5 experiments 1-8.

    The half-cycle duration is PE + BD + HP + BD-vac + LP.  Feed starts after
    PE and lasts the Table 5 feed duration.  Any remaining interval is the
    closed, heat-exchanging idle of the other column.  The second half-cycle
    swaps columns exactly as Fig. 3's caption states.
    """
    if type(experiment) is not int or experiment not in range(1, 9):
        raise ValueError("Only source-screen experiments 1-8 have complete nominal timing")
    source = json.loads(SOURCE.read_text(encoding="utf-8-sig"))
    table = source["experiment_table_5_pdf_text_screen"]
    if table["status"] != "rows_1_to_8_visually_checked_row_9_printed_conflict_preserved":
        raise ValueError("Table 5 rows 1-8 have not been source-qualified for timing")
    common = source["common_operating_table_4"]
    row = next(item for item in table["rows_1_to_8"]
               if item["experiment"] == experiment)
    pe = float(common["pressure_equalization_s"])
    bd = float(common["blowdown_s"])
    hp = float(row["heavy_purge_s"])
    vacuum = float(common["vacuum_blowdown_s"])
    lp = float(common["light_purge_s"])
    feed = float(row["feed_s"])
    if not all(math.isfinite(value) and value > 0 for value in (pe, bd, hp, vacuum, lp, feed)):
        raise ValueError("All reported operating durations must be finite and positive")
    half = pe + bd + hp + vacuum + lp
    feed_end = pe + feed
    if not (pe < feed_end < half):
        raise ValueError("The source half-cycle needs a positive feed and idle interval")
    idle = half - feed_end
    charging_boundaries = (0.0, pe, feed_end, half)
    charging_names = ("PE-Press", "Feed", "Idle")
    discharging_boundaries = (0.0, pe, pe + bd, pe + bd + hp,
                              pe + bd + hp + vacuum, half)
    discharging_names = ("PE-BD", "BD", "HP", "BD-vac", "LP")
    half_events = sorted(set(charging_boundaries + discharging_boundaries))
    slots = []
    for half_index in (0, 1):
        offset = half_index * half
        for i in range(len(half_events) - 1):
            local_start, local_end = half_events[i:i + 2]
            midpoint = (local_start + local_end) / 2.0
            charge = _operation_at(midpoint, charging_boundaries, charging_names)
            discharge = _operation_at(midpoint, discharging_boundaries, discharging_names)
            actions = [charge, discharge] if half_index == 0 else [discharge, charge]
            pairs = [[1, 0]] if half_index == 0 and local_start == 0 else (
                [[0, 1]] if half_index == 1 and local_start == 0 else []
            )
            slots.append({
                "event_interval_index": len(slots),
                "half_cycle_index": half_index,
                "start_s": offset + local_start,
                "end_s": offset + local_end,
                "duration_s": local_end - local_start,
                "actions": actions,
                "equalization_donor_receiver_pairs": pairs,
            })
    per_bed = [_collapsed_actions(slots, bed) for bed in range(2)]
    expected_sequences = [
        list(NAMED_SEQUENCE),
        list(NAMED_SEQUENCE[3:] + NAMED_SEQUENCE[:3]),
    ]
    if per_bed != expected_sequences:
        raise ArithmeticError("Merged event intervals changed the source bed sequence")
    if experiment == 1:
        checked = source["schedule_text_reconstruction"]
        if not (half == checked["experiment_1_half_cycle_s"]
                and feed_end == checked["experiment_1_feed_end_s"]
                and idle == checked["experiment_1_idle_duration_s"]):
            raise ArithmeticError("Experiment 1 diverges from the visually checked Figure 3 contract")
    return {
        "schema": "dynamic-psa-streb-two-bed-event-timeline/v1",
        "status": ("figure_3_experiment_1_timing_checked_not_engine_qualified" if experiment == 1
                   else "table_5_timing_visually_checked_not_engine_qualified"),
        "experiment": experiment,
        "number_of_beds": 2,
        "named_operations_per_bed": 8,
        "global_event_interval_count": len(slots),
        "half_cycle_time_s": half,
        "full_cycle_time_s": 2.0 * half,
        "feed_duration_s": feed,
        "idle_duration_s": idle,
        "press_ads_split_s": None,
        "bed_named_sequences": per_bed,
        "event_intervals": slots,
        "source_figure": "https://link.springer.com/article/10.1007/s10450-021-00308-w/figures/3",
        "engine_input_allowed": False,
        "release_flags": {
            "source_validation_released": False,
            "design_eligible": False,
            "dwsim_application_allowed": False,
            "tea_lca_allowed": False,
            "training_eligible": False,
        },
    }
