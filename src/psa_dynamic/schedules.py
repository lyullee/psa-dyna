"""Schedule plugins. A slot advances every bed on the same physical clock."""

import numpy as np

from .streb_2021_timeline import build_streb_two_bed_timeline


def validate_schedule_contract(schedule):
    """Validate the topology/count contract used by the dynamic runner.

    ``number_of_beds`` is the number of simultaneous action rows (physical
    beds).  ``number_of_cycle_steps`` is the number of global time slots in
    one cycle.  They are intentionally independent quantities; neither is
    inferred from the other.  This validator is also used for hand-authored
    diagnostic schedules, so it checks structure rather than requiring a
    particular source plugin.
    """
    if not isinstance(schedule, dict):
        raise TypeError("schedule must be a mapping")
    try:
        number_of_beds = int(schedule["number_of_beds"])
        number_of_cycle_steps = int(schedule["number_of_cycle_steps"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("schedule must declare integer number_of_beds and number_of_cycle_steps") from exc
    if number_of_beds < 1 or number_of_cycle_steps < 1:
        raise ValueError("number_of_beds and number_of_cycle_steps must be positive")
    if isinstance(schedule["number_of_beds"], bool) or number_of_beds != schedule["number_of_beds"]:
        raise ValueError("number_of_beds must be an integer, not a coercible value")
    if isinstance(schedule["number_of_cycle_steps"], bool) or number_of_cycle_steps != schedule["number_of_cycle_steps"]:
        raise ValueError("number_of_cycle_steps must be an integer, not a coercible value")
    slots = schedule.get("slots")
    if not isinstance(slots, list) or len(slots) != number_of_cycle_steps:
        raise ValueError(
            "cycle-step count mismatch: number_of_cycle_steps must equal len(slots)"
        )
    durations = []
    for index, slot in enumerate(slots):
        if not isinstance(slot, dict):
            raise ValueError(f"slot {index} must be a mapping")
        if slot.get("slot_index") != index:
            raise ValueError(f"slot {index} has inconsistent slot_index")
        actions = slot.get("actions")
        if not isinstance(actions, list) or len(actions) != number_of_beds:
            raise ValueError(
                f"bed count mismatch at slot {index}: actions must contain exactly "
                f"number_of_beds={number_of_beds} entries"
            )
        try:
            duration = float(slot["duration_s"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"slot {index} must declare numeric duration_s") from exc
        if not np.isfinite(duration) or duration <= 0:
            raise ValueError(f"slot {index} duration_s must be finite and positive")
        durations.append(duration)
    try:
        cycle_time_s = float(schedule["cycle_time_s"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("schedule must declare numeric cycle_time_s") from exc
    if not np.isfinite(cycle_time_s) or cycle_time_s <= 0:
        raise ValueError("cycle_time_s must be finite and positive")
    if not np.isclose(cycle_time_s, sum(durations), rtol=0.0, atol=1e-9):
        raise ValueError("cycle_time_s must equal the sum of slot durations")
    # Optional provenance markers make the independence explicit in serialized
    # schedules while remaining compatible with older hand-authored schedules.
    if "bed_count_definition" in schedule and schedule["bed_count_definition"] != "action_rows":
        raise ValueError("bed_count_definition must be 'action_rows'")
    if "cycle_step_count_definition" in schedule and schedule["cycle_step_count_definition"] != "slot_count":
        raise ValueError("cycle_step_count_definition must be 'slot_count'")
    contract = {
        "number_of_beds": number_of_beds,
        "number_of_cycle_steps": number_of_cycle_steps,
        "bed_count_definition": "action_rows",
        "cycle_step_count_definition": "slot_count",
        "counts_are_independent": True,
    }
    # Some publications call simultaneously executed sub-operations "steps".
    # Preserve that source terminology separately from the integrator's global
    # clock slots so a two-bed/six-operation source cycle is never mislabeled
    # as either six beds or six serial time intervals.
    source_operations = schedule.get("source_named_operations_per_bed")
    if source_operations is not None:
        if (isinstance(source_operations, bool)
                or int(source_operations) != source_operations
                or int(source_operations) < 1):
            raise ValueError("source_named_operations_per_bed must be a positive integer")
        contract["source_named_operations_per_bed"] = int(source_operations)
        contract["source_operation_count_definition"] = "distinct_named_operations_executed_by_each_bed"
    return contract


def build_schedule(config):
    name = config["name"]
    if name == "skarstrom_4step":
        # Two beds alternate feed pressurization/adsorption and blowdown/purge.
        rows = [["FP", "AD", "BD", "P"], ["BD", "P", "FP", "AD"]]
        durations = [config["pressurization_s"], config["adsorption_s"]]*2
        used_timing_inputs = ("pressurization_s", "adsorption_s")
        pairs = {}
        evidence = "Classical Skarstrom sequence; declared engineering timings."
    elif name == "equalized_8step":
        # Four beds, two consecutive adsorption slots, two PE donations,
        # blowdown, purge and two PE receptions. Offsets 0/1/4/5 make the two
        # donations go to different receivers. This is an engineering plugin,
        # not an exact copy of Krótki's asynchronous valve schedule.
        cycle = ["AD", "AD", "ED1", "ED2", "BD", "P", "ER2", "ER1"]
        rows = [[cycle[(slot-offset) % 8] for slot in range(8)] for offset in (0, 1, 4, 5)]
        durations = [config["slot_s"]]*8
        used_timing_inputs = ("slot_s",)
        pairs = {0: [(3, 1)], 2: [(0, 2)], 3: [(1, 2), (0, 3)],
                 4: [(1, 3)], 6: [(2, 0)], 7: [(3, 0), (2, 1)]}
        evidence = "Engineering four-bed eight-slot/two-equalization schedule; requires independent validation."
    elif name == "linde_type_8step_source":
        # Krótki et al. (Energies 19, 2753, 2026), Table 6.  The source has
        # two different inter-bed transfers: pressure equalization and a
        # product-purge transfer.  PPD/PPR are deliberately separate from
        # the historical P action, which combines counter-current tail flow
        # with an automatic purge draw.  This prevents silently changing the
        # source boundary semantics.
        rows = [
            ["AD", "IDLE", "ED2", "PPD", "BD", "PPR", "ER1", "FP"],
            ["ER1", "FP", "AD", "IDLE", "ED2", "PPD", "BD", "PPR"],
            ["BD", "PPR", "ER1", "FP", "AD", "IDLE", "ED2", "PPD"],
            ["ED2", "PPD", "BD", "PPR", "ER1", "FP", "IDLE", "AD"],
        ]
        # The article publishes installation-specific ranges for the eight
        # operations, but not a machine-readable valve-time history. Keep the
        # equal-slot default and allow an explicit eight-value sensitivity
        # vector without treating it as source timing recovery.
        if "slot_durations_s" in config:
            durations = [float(value) for value in config["slot_durations_s"]]
            if len(durations) != 8:
                raise ValueError("slot_durations_s must contain exactly eight values")
            used_timing_inputs = ("slot_durations_s",)
        else:
            durations = [config["slot_s"]] * 8
            used_timing_inputs = ("slot_s",)
        pairs = {
            0: [(3, 1)],
            2: [(0, 2)],
            4: [(1, 3)],
            6: [(2, 0)],
        }
        product_purge_pairs = {
            1: [(3, 2)],
            3: [(0, 3)],
            5: [(1, 0)],
            7: [(2, 0)],
        }
        evidence = "Krótki et al. Table 6 source topology; product-purge valve law and material data remain unbound."
    elif name == "wang_type_i_12step":
        # Transcription: Wang et al. CEJ 472 (2023) 144920, Fig. 2/S8.
        rows = [
            ["FP", "AD", "ED1", "ED2", "IDLE", "BD", "BD", "P", "IDLE", "ER2", "IDLE", "ER1"],
            ["ER2", "IDLE", "ER1", "FP", "AD", "ED1", "ED2", "IDLE", "BD", "BD", "P", "IDLE"],
            ["BD", "P", "IDLE", "ER2", "IDLE", "ER1", "FP", "AD", "ED1", "ED2", "IDLE", "BD"],
            ["ED2", "IDLE", "BD", "BD", "P", "IDLE", "ER2", "IDLE", "ER1", "FP", "AD", "ED1"],
        ]
        durations = [config["slot_s"]]*12
        used_timing_inputs = ("slot_s",)
        pairs = {0:[(3,1)], 2:[(0,1)], 3:[(0,2)], 5:[(1,2)], 6:[(1,3)], 8:[(2,3)], 9:[(2,0)], 11:[(3,0)]}
        evidence = "Published topology; equal-slot timing is a pressure-figure inference, not a reported valve log."
    elif name == "welter_reduced_3step_protocol":
        # Welter et al. (2026), DOI 10.1002/cite.70127: one physical bed
        # executes adsorption, depressurization and evacuation for equal
        # durations. EVC is intentionally not represented by the historical
        # P action: the experiment has no simultaneous adsorption donor or
        # light-product purge. This plugin closes only the reported topology;
        # source-equivalent timing remains event-driven and the total-pressure
        # boundary still requires a hash-pinned source trace.
        rows = [["AD", "BD", "EVC"]]
        durations = [config["step_s"]] * 3
        used_timing_inputs = ("step_s",)
        pairs = {}
        evidence = (
            "Welter et al. 2026 one-bed AD-BD-EVC protocol topology; "
            "fixed step_s is a numerical placeholder, not source-equivalent event timing."
        )
    elif name == "streb_two_bed_8operation_screen":
        # Figure 3 has eight named operations per bed, but their transitions
        # create more than eight global event intervals.  This adapter is a
        # topology diagnostic and requires the dedicated Network boundary
        # mode; it is not the paper's pressure trajectory or RAST replay.
        timeline = build_streb_two_bed_timeline(config["experiment"])
        rows = [[slot["actions"][bed] for slot in timeline["event_intervals"]]
                for bed in range(timeline["number_of_beds"])]
        durations = [slot["duration_s"] for slot in timeline["event_intervals"]]
        pairs = {index: slot["equalization_donor_receiver_pairs"]
                 for index, slot in enumerate(timeline["event_intervals"])
                 if slot["equalization_donor_receiver_pairs"]}
        used_timing_inputs = ("experiment",)
        evidence = "Streb-Mazzotti Part 2 Figure 3/Table 5 timing; diagnostic boundary topology only."
    elif name == "guderian_two_bed_6step_source":
        # Marcinek et al. (2021), Fig. 8/Table 5, and Marcinek's dissertation
        # Table 4.3-2 distinguish six *named operations per bed* from four
        # simultaneous global clock intervals.  During each 59-s process
        # interval the low-pressure equalization receiver performs feed/
        # product pressurisation followed by production, while the original
        # high-pressure donor performs blowdown with receiver purge.  Keeping
        # that phase relation explicit is essential: reversing AP/BR after
        # equalization produces the wrong rising/falling pressure traces even
        # though the four interval durations still look valid.
        # The transition within each composite operation is pressure/direction
        # driven and must not be fabricated as extra fixed-duration slots.
        half_cycle = float(config["half_cycle_s"])
        equalization = float(config["equalization_s"])
        process = half_cycle - equalization
        if process <= 0.0:
            raise ValueError("Guderian half_cycle_s must exceed equalization_s")
        rows = [
            ["EQ-D", "BR", "EQ-R", "AP"],
            ["EQ-R", "AP", "EQ-D", "BR"],
        ]
        durations = [equalization, process, equalization, process]
        used_timing_inputs = ("half_cycle_s", "equalization_s")
        # Top/top and bottom/bottom connections operate simultaneously.  The
        # legacy equalization_pairs field is intentionally empty because it
        # can represent only a product-end transfer.  The source adapter must
        # consume the endpoint-explicit connection table below.
        pairs = {}
        guderian_connections = {
            0: [
                {"donor": 0, "receiver": 1, "donor_endpoint": "top", "receiver_endpoint": "top"},
                {"donor": 0, "receiver": 1, "donor_endpoint": "bottom", "receiver_endpoint": "bottom"},
            ],
            2: [
                {"donor": 1, "receiver": 0, "donor_endpoint": "top", "receiver_endpoint": "top"},
                {"donor": 1, "receiver": 0, "donor_endpoint": "bottom", "receiver_endpoint": "bottom"},
            ],
        }
        evidence = (
            "Marcinek et al. 2021 Fig. 8/Table 5 plus Marcinek dissertation "
            "Table 4.3-2: 1-s equalization and 59-s simultaneous production/"
            "blowdown/purge at the 60-s reference half-cycle."
        )
    elif name == "breakthrough":
        rows = [["AD"]]
        durations = [config["duration_s"]]
        used_timing_inputs = ("duration_s",)
        pairs = {}
        evidence = "Single continuous breakthrough experiment; CSS is not applicable."
    elif name == "closed":
        rows = [["IDLE"]]
        durations = [config["duration_s"]]
        used_timing_inputs = ("duration_s",)
        pairs = {}
        evidence = "Verification-only isolated bed."
    else:
        raise ValueError(f"Unknown schedule: {name}")
    if any(d <= 0 for d in durations):
        raise ValueError("Slot durations must be positive")
    # Historical case files may contain physical timing fields which a plugin
    # does not use.  Expose them in every built schedule so downstream audits
    # cannot mistake an ignored value for an executed valve duration.
    timing_keys = {
        "slot_s", "slot_durations_s", "pressurization_s", "adsorption_s",
        "step_s", "half_cycle_s", "equalization_s",
    }
    if name in {"breakthrough", "closed"}:
        timing_keys.add("duration_s")
    ignored_timing_inputs = {
        key: config[key] for key in sorted(timing_keys - set(used_timing_inputs))
        if key in config
    }
    # A cyclic case may carry a legacy ``duration_s`` field inherited from a
    # breakthrough case.  Do not interpret that value as a cycle duration
    # without source timing evidence.  Only the explicit
    # ``declared_cycle_time_s`` field is treated as a reported cycle period.
    # One-slot breakthrough/closed cases still use ``duration_s`` physically.
    cycle_time_s = float(sum(durations))
    cyclic_schedule = name in {
        "skarstrom_4step", "equalized_8step", "linde_type_8step_source",
        "wang_type_i_12step", "streb_two_bed_8operation_screen",
        "welter_reduced_3step_protocol", "guderian_two_bed_6step_source",
    }
    legacy_duration_s = config.get("duration_s") if cyclic_schedule else None
    reported_duration_s = config.get("reported_breakthrough_duration_s") if cyclic_schedule else None
    declared_cycle_time_s = config.get("declared_cycle_time_s")
    duration_difference_s = (
        None if declared_cycle_time_s is None
        else float(cycle_time_s - float(declared_cycle_time_s))
    )
    cycle_time_consistent = (
        True if duration_difference_s is None
        else abs(duration_difference_s) <= 1e-9
    )
    if declared_cycle_time_s is not None:
        duration_role = "declared_cycle_time"
    elif legacy_duration_s is not None:
        duration_role = "legacy_duration_not_assumed_cycle"
    elif reported_duration_s is not None:
        duration_role = "reported_noncycle_duration"
    elif name in {"breakthrough", "closed"}:
        duration_role = "single_slot_experiment_duration"
    else:
        duration_role = "inferred_from_slot_durations"
    slots = [{"slot_index":i, "duration_s":float(d), "actions":[row[i] for row in rows],
              "equalization_pairs":pairs.get(i, []),
              "product_purge_pairs":product_purge_pairs.get(i, []) if name == "linde_type_8step_source" else [],
              "equalization_connections":guderian_connections.get(i, []) if name == "guderian_two_bed_6step_source" else []}
             for i, d in enumerate(durations)]
    for slot in slots:
        involved = []
        for donor, receiver in slot["equalization_pairs"]:
            if name == "streb_two_bed_8operation_screen":
                pair_matches = (slot["actions"][donor] == "PE-BD"
                                and slot["actions"][receiver] == "PE-Press")
            else:
                pair_matches = (slot["actions"][donor].startswith("ED")
                                and slot["actions"][receiver].startswith("ER"))
            if not pair_matches:
                raise ValueError("Equalization connection does not match the valve table")
            involved.extend([donor, receiver])
        if len(involved) != len(set(involved)):
            raise ValueError("A bed cannot be in two simultaneous equalizations")
        if name == "streb_two_bed_8operation_screen":
            donors = sum(a == "PE-BD" for a in slot["actions"])
            receivers = sum(a == "PE-Press" for a in slot["actions"])
        else:
            donors = sum(a.startswith("ED") for a in slot["actions"])
            receivers = sum(a.startswith("ER") for a in slot["actions"])
        if donors != len(slot["equalization_pairs"]) or receivers != donors:
            raise ValueError("Unconnected equalization operation")
        source_connections = slot.get("equalization_connections", [])
        if source_connections:
            source_donors = [index for index, action in enumerate(slot["actions"]) if action == "EQ-D"]
            source_receivers = [index for index, action in enumerate(slot["actions"]) if action == "EQ-R"]
            if len(source_donors) != 1 or len(source_receivers) != 1 or len(source_connections) != 2:
                raise ValueError("Guderian equalization requires one donor, one receiver and two endpoint connections")
            endpoints = set()
            for connection in source_connections:
                if (connection.get("donor") != source_donors[0]
                        or connection.get("receiver") != source_receivers[0]):
                    raise ValueError("Guderian equalization connection does not match the action rows")
                endpoint_pair = (connection.get("donor_endpoint"), connection.get("receiver_endpoint"))
                if endpoint_pair not in {("top", "top"), ("bottom", "bottom")}:
                    raise ValueError("Guderian equalization must connect like endpoints")
                endpoints.add(endpoint_pair)
            if endpoints != {("top", "top"), ("bottom", "bottom")}:
                raise ValueError("Guderian equalization requires simultaneous top and bottom connections")
        product_pairs = slot.get("product_purge_pairs", [])
        for donor, receiver in product_pairs:
            # In Table 6, the final 3 -> 1 product-purge transfer enters the
            # receiving bed while it is compressed with product rather than a
            # separate PPR-labelled state.
            if slot["actions"][donor] != "PPD" or slot["actions"][receiver] not in {"PPR", "FP"}:
                raise ValueError("Product-purge connection does not match the source valve table")
        product_donors = sum(a == "PPD" for a in slot["actions"])
        product_receivers = len({receiver for _, receiver in product_pairs})
        if product_donors != len(product_pairs) or product_receivers != product_donors:
            raise ValueError("Unconnected product-purge operation")
        if "P" in slot["actions"] and "AD" not in slot["actions"]:
            raise ValueError("This plugin needs simultaneous adsorption to supply purge")
    if name == "guderian_two_bed_6step_source":
        # A pressure-equalization donor is the bed leaving production and must
        # enter blowdown/regeneration; the receiver is the regenerated bed and
        # must enter adsorption/production.  Validate the phase relation, not
        # only the number and duration of the four global intervals.
        for equalization_index, process_index in ((0, 1), (2, 3)):
            connections = slots[equalization_index]["equalization_connections"]
            donor = connections[0]["donor"]
            receiver = connections[0]["receiver"]
            following = slots[process_index]["actions"]
            if following[donor] != "BR" or following[receiver] != "AP":
                raise ValueError(
                    "Guderian equalization donor must enter BR and receiver must enter AP"
                )
    return {"name":name, "number_of_beds":len(rows), "number_of_cycle_steps":len(slots),
            "bed_count_definition":"action_rows",
            "cycle_step_count_definition":"slot_count",
            "counts_are_independent":True,
            "source_named_operations_per_bed":({
                "streb_two_bed_8operation_screen": 8,
                "guderian_two_bed_6step_source": 6,
            }.get(name)),
            "source_experiment":config["experiment"] if name == "streb_two_bed_8operation_screen" else None,
            "source_equivalent":False if name in {
                "streb_two_bed_8operation_screen", "welter_reduced_3step_protocol",
                "guderian_two_bed_6step_source",
            } else None,
            "used_timing_inputs":{key:config[key] for key in used_timing_inputs},
            "ignored_timing_inputs":ignored_timing_inputs,
            "timing_inputs_unambiguous":not ignored_timing_inputs,
            "cycle_time_s":cycle_time_s, "declared_cycle_time_s":(
                None if declared_cycle_time_s is None else float(declared_cycle_time_s)),
            "legacy_duration_s":(
                None if legacy_duration_s is None else float(legacy_duration_s)),
            "reported_breakthrough_duration_s":(
                None if reported_duration_s is None else float(reported_duration_s)),
            "duration_role":duration_role,
            "cycle_time_difference_s":duration_difference_s,
            "cycle_time_consistent":cycle_time_consistent,
            "slots":slots, "evidence":evidence}
