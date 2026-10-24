"""Independent, restartable case execution for local and PBS-array workers."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import time
import zipfile
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp,quad_vec

from .core import Network
from .schedules import build_schedule
from .numerics import BoundedJacobian,ComplexStepJacobian,dense_gauss5_with_gauss3_error
from .qualification import case_screen, update_trajectory_extrema, trajectory_domain_screen


# Only modules that can change the numerical state transition, schedule
# semantics, or qualification calculations belong in the restart identity.
# Audit/report scripts are intentionally excluded: editing a register or a
# literature-summary script must not invalidate a long-running numerical
# checkpoint.  The full source snapshot is still archived separately by
# ``snapshot_code`` for provenance.
ENGINE_SOURCE_FILES = (
    "core.py",
    "guderian_2021_n2.py",
    "numerics.py",
    "qualification.py",
    "runner.py",
    "schedules.py",
    "streb_2021_pressure.py",
    "streb_2021_timeline.py",
    "streb_2021_wilson.py",
    "vectorized.py",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_hash():
    checksum = hashlib.sha256()
    package_root = Path(__file__).parent
    for name in ENGINE_SOURCE_FILES:
        path = package_root / name
        if not path.is_file():
            raise FileNotFoundError(f"Required engine source is missing: {path}")
        checksum.update(path.name.encode()); checksum.update(path.read_bytes())
    return checksum.hexdigest()


def snapshot_code(output):
    path=output/"source_code.zip"
    if not path.exists():
        temporary=output/"source_code.zip.tmp"
        with zipfile.ZipFile(temporary,"w",compression=zipfile.ZIP_DEFLATED) as archive:
            for source in sorted(Path(__file__).parent.glob("*.py")):
                info=zipfile.ZipInfo("psa_dynamic/"+source.name,date_time=(1980,1,1,0,0,0))
                info.compress_type=zipfile.ZIP_DEFLATED
                archive.writestr(info,source.read_bytes())
        os.replace(temporary,path)
    return digest(path)


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.flush(); os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save_checkpoint(directory, metadata, state):
    # One compressed array file per case, not thousands of timestep files.
    # Alternate generations. If a job dies before JSON commit, the preceding
    # JSON still points to an untouched, complete array generation.
    generation = 1-int(metadata.get("state_generation",0))
    array_path = directory/f"checkpoint.state.{generation}.npz"
    temporary = array_path.with_suffix(".npz.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, state=state)
        handle.flush(); os.fsync(handle.fileno())
    os.replace(temporary, array_path)
    metadata["state_generation"] = generation
    metadata["state_file"] = array_path.name
    metadata["state_sha256"] = digest(array_path)
    atomic_json(directory/"checkpoint.json", metadata)


def archive_failure(directory):
    """Preserve a previous caught failure before a hash-verified new attempt."""
    failure=directory/"failure.json"
    if not failure.exists():
        return
    history=directory/"failure_history"
    history.mkdir(exist_ok=True)
    index=1
    while (history/f"failure_{index:06d}.json").exists():
        index+=1
    atomic_json(history/f"failure_{index:06d}.json",json.loads(failure.read_text()))
    failure.unlink()


def scaled_state_change(network, before, after):
    ns = network.ns
    left = before[:network.bed_process_size].reshape(network.nb, network.bed.n, network.bed.nv)
    right = after[:network.bed_process_size].reshape(left.shape)
    floor = np.r_[np.full(ns, 1e-8), np.full(ns, 1e-4), 1., 1.]
    scale = np.maximum(np.maximum(abs(left), abs(right)), floor)
    scaled=abs(right-left)/scale
    worst=np.unravel_index(np.argmax(scaled),scaled.shape)
    names=["gas_"+s["name"] for s in network.species]+[
        ("pellet_pressure_" if network.bed.mode=="partial_pressure_ldf" else "loading_")+s["name"] for s in network.species]+["bed_temperature","wall_temperature"]
    result = {"state_scaled_max":float(scaled[worst]),"state_worst_bed_index":int(worst[0]),
              "state_worst_cell_index":int(worst[1]),"state_worst_variable":names[worst[2]],
              "state_worst_before":float(left[worst]),"state_worst_after":float(right[worst])}
    pressures = []
    temperatures = []
    for b in range(network.nb):
        lb = network.bed.fields(before[b*network.bed.size:(b+1)*network.bed.size])
        rb = network.bed.fields(after[b*network.bed.size:(b+1)*network.bed.size])
        pressures.append(np.max(abs(lb[4]-rb[4])))
        temperatures.append(np.max(abs(lb[2]-rb[2])))
    result["pressure_max_abs_bar"] = float(max(pressures))
    result["temperature_max_abs_k"] = float(max(temperatures))
    if network.receiver is not None:
        receiver_left = before[network.receiver_offset:network.receiver_offset+network.receiver_size]
        receiver_right = after[network.receiver_offset:network.receiver_offset+network.receiver_size]
        receiver_floor = np.r_[np.full(ns, 1e-8), 1.0]
        receiver_scale = np.maximum(np.maximum(abs(receiver_left), abs(receiver_right)), receiver_floor)
        receiver_scaled = abs(receiver_right-receiver_left)/receiver_scale
        receiver_worst = int(np.argmax(receiver_scaled))
        if receiver_scaled[receiver_worst] > result["state_scaled_max"]:
            result.update({
                "state_scaled_max": float(receiver_scaled[receiver_worst]),
                "state_worst_bed_index": -1,
                "state_worst_cell_index": -1,
                "state_worst_variable": (
                    f"receiver_inventory_{network.species[receiver_worst]['name']}"
                    if receiver_worst < ns else "receiver_temperature"
                ),
                "state_worst_before": float(receiver_left[receiver_worst]),
                "state_worst_after": float(receiver_right[receiver_worst]),
            })
        left_fields = network.receiver.fields(receiver_left)
        right_fields = network.receiver.fields(receiver_right)
        result["receiver_pressure_abs_change_bar"] = float(abs(left_fields[2]-right_fields[2]))
        result["receiver_temperature_abs_change_k"] = float(abs(left_fields[1]-right_fields[1]))
        result["pressure_max_abs_bar"] = max(
            result["pressure_max_abs_bar"], result["receiver_pressure_abs_change_bar"]
        )
        result["temperature_max_abs_k"] = max(
            result["temperature_max_abs_k"], result["receiver_temperature_abs_change_k"]
        )
    return result


def run_case(case_path, output, resume=True, cycle_budget=None, step_budget=None):
    """Run one case. A paused/failed/unconverged case is never training truth.

    The cluster invokes this same function once per immutable input JSON.
    Interrupted jobs resume at the last completed physical slot. A hash change
    rejects the checkpoint instead of silently mixing model versions.
    """
    case_path, output = Path(case_path), Path(output)
    case = json.loads(case_path.read_text(encoding="utf-8-sig"))
    network = Network(case, build_schedule(case["schedule"]))
    output.mkdir(parents=True, exist_ok=True)
    lock_path = output/"RUNNING.lock"
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f"Case has a running lock: {lock_path}. Inspect the owning job before removing it.")
    os.write(fd, f"{platform.node()} pid={os.getpid()}\n".encode()); os.close(fd)
    hashes = {"input_sha256":digest(case_path), "code_sha256":code_hash()}
    started = time.perf_counter()
    metadata = {}
    try:
        result_path = output/"result.json"
        if result_path.exists() and resume:
            cached = json.loads(result_path.read_text())
            if cached.get("status") in {"css_converged", "experiment_completed"}:
                if any(cached.get(k) != v for k,v in hashes.items()):
                    raise ValueError("Cached result hash differs from requested code/input")
                return cached
        checkpoint = output/"checkpoint.json"
        if checkpoint.exists() and resume:
            metadata = json.loads(checkpoint.read_text())
            if any(metadata.get(k) != v for k,v in hashes.items()):
                raise ValueError("Checkpoint code/input hashes differ; use a fresh output directory")
            if digest(output/metadata["state_file"]) != metadata["state_sha256"]:
                raise ValueError("Checkpoint array hash mismatch")
            with np.load(output/metadata["state_file"], allow_pickle=False) as archive:
                state = archive["state"]
            metadata["resume_count"] += 1
        else:
            state = network.initial()
            if "initial_process_state" in case:
                seed=np.asarray(case["initial_process_state"],dtype=float)
                if seed.shape!=(network.process_size,) or not np.isfinite(seed).all():
                    raise ValueError("Initial process-state seed has incompatible dimensions")
                # A warm start is an initial condition, never a CSS assertion.
                # Initial bed inventory is recomputed and still enters closure.
                state[:network.process_size]=seed
            inv, energy = network.inventory_energy(state)
            metadata = {**hashes, "schema":"dynamic-psa-checkpoint/v1", "cycle_index":0,"slot_index":0,
                        "initial_inventory":inv.tolist(),"initial_energy_j":energy,"resume_count":0,
                        "cycle_start_state":state.tolist(),"cycle_rows":[],"step_rows":[],"trace_rows":[],"css_streak":0,
                        "wall_seconds_total":0.,"trajectory_extrema":{}}
        initial_wall_seconds = metadata["wall_seconds_total"]
        archive_failure(output)
        atomic_json(output/"input.json", case)
        atomic_json(output/"schedule.json", network.schedule)
        input_screen = case_screen(case)
        atomic_json(output/"engineering_input_screen.json",input_screen)
        source_zip_hash=snapshot_code(output)
        solver = case["solver"]
        atol = np.tile(np.r_[np.full(network.ns, 1e-12),np.full(network.ns,2e-10),3e-7,3e-7], network.nb*network.bed.n)
        if network.receiver is not None:
            atol = np.r_[atol, np.full(network.ns, 1e-12), 3e-7]
        atol = np.r_[atol, np.full(3*network.ns,1e-11),np.full(5,1e-7)]
        sparsity = network.sparsity()
        cycles_done = 0; slots_done = 0; status = "maximum_cycles_reached"
        noncyclic = network.schedule["name"] in {"breakthrough", "closed"}
        maximum_cycles = 1 if noncyclic else solver["maximum_cycles"]
        while metadata["cycle_index"] < maximum_cycles:
            cycle = metadata["cycle_index"]
            for slot_index in range(metadata["slot_index"], network.schedule["number_of_cycle_steps"]):
                slot = network.schedule["slots"][slot_index]
                stamp = time.perf_counter()
                record_trace = solver.get("record_trace", False)
                record_boundary_diagnostics = solver.get("record_boundary_diagnostics", False)
                ledger_mode=solver.get("ledger_integration","quadrature")
                if ledger_mode not in {"quadrature","augmented","gauss5_trial"}:
                    raise ValueError("Unknown stream-ledger integration mode")
                coupled=ledger_mode=="augmented"
                before_ledger=state[network.process_size:].copy()
                slot_pattern=network.sparsity(slot)
                if coupled:
                    function=lambda t,x:network.rhs(t,x,slot)
                    initial=state
                    active_atol=atol
                    active_pattern=slot_pattern
                else:
                    function=lambda t,x:network.rhs(t,np.r_[x,np.zeros(network.nledger)],slot)[:network.process_size]
                    initial=state[:network.process_size]
                    active_atol=atol[:network.process_size]
                    active_pattern=slot_pattern[:network.process_size,:network.process_size]
                raw_function=function
                progress={"last_report_wall":stamp,"evaluations":0}
                def guarded_function(t,x):
                    now=time.perf_counter()
                    progress["evaluations"]+=1
                    elapsed=now-stamp
                    if now-progress["last_report_wall"]>30:
                        atomic_json(output/"active_step.json",{"status":"running","cycle":cycle+1,"slot":slot_index+1,
                            "physical_time_in_slot_s":float(t),"slot_duration_s":slot["duration_s"],
                            "wall_seconds":elapsed,"rhs_evaluations":progress["evaluations"]})
                        progress["last_report_wall"]=now
                    if elapsed>solver.get("maximum_step_wall_s",600.):
                        raise TimeoutError(f"Step wall budget: cycle={cycle+1}, slot={slot_index+1}, physical t={t:.8g} s, wall={elapsed:.1f} s; no performance/CSS claim")
                    return raw_function(t,x)
                function=guarded_function
                # Complex-step differentiates the current active physical
                # branch without subtracting nearly equal RHS evaluations.
                # Bounded and SciPy FD remain independent cross-checks.
                jacobian_mode=solver.get("jacobian_method","complex_step")
                jacobian_options={"jac_sparsity":active_pattern}
                if jacobian_mode in {"bounded_fd","complex_step"}:
                    gas_scale=case["high_pressure_bar"]*1e5*network.bed.vgas/(8.314462618*case["bed"]["ambient_k"])
                    ads_scale=case["high_pressure_bar"] if network.bed.mode=="partial_pressure_ldf" else 1.
                    one_bed=np.concatenate([np.r_[np.full(network.ns,gas),np.full(network.ns,ads_scale),300.,300.] for gas in gas_scale])
                    scales=np.tile(one_bed,network.nb)
                    if network.receiver is not None:
                        receiver_inventory = np.asarray(initial)[
                            network.receiver_offset:network.receiver_offset+network.ns
                        ]
                        scales = np.r_[scales, np.maximum(np.abs(receiver_inventory), 1e-6), 300.0]
                    if coupled:
                        scales=np.r_[scales,np.full(network.nledger,1.)]
                    jacobian_type=ComplexStepJacobian if jacobian_mode=="complex_step" else BoundedJacobian
                    jacobian_options={"jac":jacobian_type(function,active_pattern,scales,solver.get("jacobian_relative_step"))}
                elif jacobian_mode!="scipy_fd":
                    raise ValueError("Unknown Jacobian method")
                solution = solve_ivp(function, (0.,slot["duration_s"]), initial,
                                     method="BDF",rtol=solver["rtol"],atol=active_atol,
                                     max_step=solver["max_step_s"],
                                     dense_output=(record_trace or not coupled or record_boundary_diagnostics),
                                     **jacobian_options)
                if not solution.success:
                    raise RuntimeError(f"cycle={cycle+1}, slot={slot_index+1}: {solution.message}")
                boundary_integrals = None
                boundary_integral_error = None
                ledger_fast_fallback = False
                if record_boundary_diagnostics:
                    # Boundary diagnostics are normally endpoint-only.  A
                    # phase-aligned audit can request the same rates integrated
                    # over the accepted BDF trajectory.  This keeps purge,
                    # product, and valve semantics separate from the physical
                    # state equations while preventing a single CSS snapshot
                    # from being misread as a cycle-average flow.
                    def boundary_vector(time_value):
                        sample = np.asarray(solution.sol(time_value))
                        process_sample = sample[:network.process_size] if coupled else sample
                        diagnostics = {}
                        network.rhs(
                            time_value,
                            np.r_[process_sample, np.zeros(network.nledger)],
                            slot,
                            ledger_only=True,
                            diagnostics=diagnostics,
                        )
                        aggregate = [
                            np.sum(diagnostics["feed_flow_by_bed_mol_s"]),
                            np.sum(diagnostics["adsorption_gross_flow_by_bed_mol_s"]),
                            np.sum(diagnostics["tail_flow_by_bed_mol_s"]),
                            np.sum(diagnostics["purge_received_flow_by_bed_mol_s"]),
                            np.sum(diagnostics["purge_donated_flow_by_bed_mol_s"]),
                            (np.sum(diagnostics["product_exported_flow_by_bed_mol_s"])
                             + diagnostics.get("receiver_product_export_mol_s", 0.0)),
                            diagnostics["purge_requested_mol_s"],
                            diagnostics["purge_unmet_mol_s"],
                            np.sum([item["flow_mol_s"] for item in diagnostics["equalization_flow_by_pair_mol_s"]]),
                        ]
                        equalization_species = np.sum(
                            [item["species_flow_mol_s"] for item in diagnostics["equalization_flow_by_pair_mol_s"]],
                            axis=0,
                        ) if diagnostics["equalization_flow_by_pair_mol_s"] else np.zeros(network.ns)
                        # Optional diagnostic detail only: integrate the feed
                        # delivered to each bed on the same accepted trajectory
                        # and the internal equalization species transfer.
                        # The first nine aggregate entries remain unchanged.
                        return np.asarray(aggregate + diagnostics["feed_flow_by_bed_mol_s"]
                                          + equalization_species.tolist(), dtype=float)

                    boundary_integrals, boundary_integral_error = quad_vec(
                        boundary_vector,
                        0.0,
                        slot["duration_s"],
                        points=solution.t[1:-1],
                        epsabs=solver.get("ledger_epsabs", 1e-9),
                        epsrel=solver.get("ledger_epsrel", 1e-9),
                        limit=max(1000, 4 * len(solution.t)),
                    )
                if coupled:
                    state=solution.y[:,-1]
                    quadrature_error=None
                else:
                    # Cumulative ledgers never affect bed dynamics. Integrate
                    # them on the actual BDF dense trajectory, subdivided at
                    # all accepted solver times, with adaptive Gauss-Kronrod.
                    ledger_function = lambda t: network.rhs(
                        t, np.r_[solution.sol(t), np.zeros(network.nledger)],
                        slot, ledger_only=True)
                    if ledger_mode == "gauss5_trial":
                        increments, estimated_errors = dense_gauss5_with_gauss3_error(
                            ledger_function, solution.t)
                        allowed_errors = (solver.get("ledger_epsabs", 1e-9)
                                          + solver.get("ledger_epsrel", 1e-9) * np.abs(increments))
                        if np.all(estimated_errors <= allowed_errors):
                            quadrature_error = float(np.max(estimated_errors))
                        else:
                            ledger_fast_fallback = True
                            increments, quadrature_error = quad_vec(
                                ledger_function, 0., slot["duration_s"],
                                points=solution.t[1:-1],
                                epsabs=solver.get("ledger_epsabs", 1e-9),
                                epsrel=solver.get("ledger_epsrel", 1e-9),
                                limit=max(1000, 4 * len(solution.t)))
                    else:
                        increments,quadrature_error=quad_vec(
                            ledger_function, 0.,slot["duration_s"],points=solution.t[1:-1],
                            epsabs=solver.get("ledger_epsabs",1e-9),epsrel=solver.get("ledger_epsrel",1e-9),
                            limit=max(1000,4*len(solution.t)))
                    state=np.r_[solution.y[:,-1],before_ledger+increments]
                # Audit every accepted BDF node, not only each slot endpoint.
                accepted=solution.y[:network.bed_process_size].reshape(network.nb,network.bed.n,network.bed.nv,-1)
                update_trajectory_extrema(
                    metadata.setdefault("trajectory_extrema", {}),
                    accepted,
                    network.bed,
                    species_names=[species["name"] for species in case["species"]],
                )
                minimum = float(accepted[:,:,:2*network.ns,:].min())
                if network.receiver is not None:
                    receiver_minimum = float(solution.y[
                        network.receiver_offset:network.receiver_offset+network.ns
                    ].min())
                    minimum = min(minimum, receiver_minimum)
                if minimum < -solver.get("negative_tolerance",1e-8) or not np.isfinite(state).all():
                    raise FloatingPointError(f"Non-finite/negative physical state: min={minimum}")
                row = {"cycle":cycle+1,"slot":slot_index+1,"physical_time_s":cycle*network.schedule["cycle_time_s"]+sum(s["duration_s"] for s in network.schedule["slots"][:slot_index+1]),
                       "nfev":int(solution.nfev),"wall_seconds":time.perf_counter()-stamp,"minimum_inventory_or_pellet_pressure":minimum,
                       "njev":int(solution.njev),"nlu":int(solution.nlu),"jacobian_mode":jacobian_mode,
                       "jacobian_colors":len(jacobian_options["jac"].groups) if "jac" in jacobian_options else None,
                       "process_rhs_evaluations_including_jacobian":progress["evaluations"],
                       "ledger_mode":ledger_mode,"ledger_fast_fallback":ledger_fast_fallback,
                       "ledger_quadrature_error_estimate":None if quadrature_error is None else float(quadrature_error)}
                # Optional endpoint-only boundary audit.  This is intentionally
                # disabled for ordinary runs: it adds no physics and only
                # records scalar feed/product/purge/tail splits for a separate
                # source-semantics diagnostic campaign.
                if solver.get("record_boundary_diagnostics", False):
                    # The cumulative species ledgers are part of the solved
                    # state. Their slot-to-slot difference gives an exact
                    # phase attribution of exported product impurities, with
                    # no additional valve or adsorption assumption.
                    slot_ledger = state[network.process_size:] - before_ledger
                    row["boundary_integrated_feed_by_component_mol"] = json.dumps(
                        slot_ledger[:network.ns].tolist(), separators=(",", ":")
                    )
                    row["boundary_integrated_product_by_component_mol"] = json.dumps(
                        slot_ledger[network.ns:2 * network.ns].tolist(), separators=(",", ":")
                    )
                    row["boundary_integrated_tail_by_component_mol"] = json.dumps(
                        slot_ledger[2 * network.ns:3 * network.ns].tolist(), separators=(",", ":")
                    )
                    boundary_diagnostics = {}
                    network.rhs(0.0, np.r_[state[:network.process_size], np.zeros(network.nledger)],
                                slot, ledger_only=True, diagnostics=boundary_diagnostics)
                    for key, value in boundary_diagnostics.items():
                        row[f"boundary_{key}"] = json.dumps(value, separators=(",", ":")) if isinstance(value, (list, dict)) else value
                    if boundary_integrals is not None:
                        labels = (
                            "feed_mol",
                            "adsorption_gross_mol",
                            "tail_mol",
                            "purge_received_mol",
                            "purge_donated_mol",
                            "product_exported_mol",
                            "purge_requested_mol",
                            "purge_unmet_mol",
                            "equalization_mol",
                        )
                        for label, value in zip(labels, boundary_integrals):
                            row[f"boundary_integrated_{label}"] = float(value)
                        per_bed_feed = np.asarray(
                            boundary_integrals[len(labels):len(labels) + network.nb], dtype=float
                        )
                        if len(per_bed_feed) != network.nb or not np.isclose(
                            per_bed_feed.sum(), boundary_integrals[0], rtol=1e-7, atol=1e-7
                        ):
                            raise ValueError("Integrated per-bed and aggregate feed diagnostics disagree")
                        row["boundary_integrated_feed_by_bed_mol"] = json.dumps(
                            per_bed_feed.tolist(), separators=(",", ":")
                        )
                        equalization_species = np.asarray(
                            boundary_integrals[len(labels) + network.nb:], dtype=float
                        )
                        if (len(equalization_species) != network.ns
                                or not np.isclose(equalization_species.sum(),
                                                  boundary_integrals[8], rtol=1e-7, atol=1e-7)):
                            raise ValueError("Integrated equalization species and total disagree")
                        row["boundary_integrated_equalization_by_component_mol"] = json.dumps(
                            equalization_species.tolist(), separators=(",", ":")
                        )
                        row["boundary_integral_quad_error"] = float(boundary_integral_error)
                for bed in range(network.nb):
                    fields = network.bed.fields(state[bed*network.bed.size:(bed+1)*network.bed.size])
                    row[f"bed{bed+1}_pressure_bar"] = float(fields[4].mean())
                    row[f"bed{bed+1}_temperature_k"] = float(fields[2].mean())
                    row[f"bed{bed+1}_temperature_min_k"]=float(accepted[bed,:,-2,:].min())
                    row[f"bed{bed+1}_temperature_max_k"]=float(accepted[bed,:,-2,:].max())
                    row[f"bed{bed+1}_wall_temperature_min_k"]=float(accepted[bed,:,-1,:].min())
                    row[f"bed{bed+1}_wall_temperature_max_k"]=float(accepted[bed,:,-1,:].max())
                metadata["step_rows"].append(row)
                if record_trace:
                    # Keep only the latest cycle's time trace. Long sweeps
                    # otherwise create enormous, redundant shared-FS output.
                    if slot_index == 0:
                        metadata["trace_rows"] = []
                    interval=solver.get("trace_interval_s",2.)
                    times=np.linspace(0,slot["duration_s"],max(2,int(np.ceil(slot["duration_s"]/interval))+1))
                    for time_value in times:
                        sample=solution.sol(time_value)
                        physical=cycle*network.schedule["cycle_time_s"]+sum(s["duration_s"] for s in network.schedule["slots"][:slot_index])+time_value
                        trace={"cycle":cycle+1,"slot":slot_index+1,"time_s":physical}
                        for b in range(network.nb):
                            fields=network.bed.fields(sample[b*network.bed.size:(b+1)*network.bed.size])
                            trace[f"bed{b+1}_top_pressure_bar"]=float(fields[4][-1])
                            # A published pressure trace may have been sampled
                            # elsewhere along the bed. Retain the axial field
                            # for non-fitting sensor-position diagnostics;
                            # do not infer a source sensor location here.
                            for cell in range(network.bed.n):
                                trace[f"bed{b+1}_cell{cell+1}_pressure_bar"]=float(fields[4][cell])
                            for j,species in enumerate(case["species"]):
                                trace[f"bed{b+1}_outlet_{species['name']}"]=float(fields[5][-1,j])
                            for cell in range(network.bed.n):
                                trace[f"bed{b+1}_cell{cell+1}_temperature_k"]=float(fields[2][cell])
                        if network.receiver is not None:
                            receiver = network.receiver.fields(sample[
                                network.receiver_offset:network.receiver_offset+network.receiver_size
                            ])
                            trace["receiver_pressure_bar"] = float(receiver[2])
                            trace["receiver_temperature_k"] = float(receiver[1])
                            for index, species in enumerate(case["species"]):
                                trace[f"receiver_{species['name']}_mole_fraction"] = float(receiver[3][index])
                        metadata["trace_rows"].append(trace)
                metadata["slot_index"] = slot_index+1
                metadata["wall_seconds_total"] = initial_wall_seconds+time.perf_counter()-started
                save_checkpoint(output, metadata, state)
                slots_done += 1
                if step_budget is not None and slots_done >= step_budget:
                    status = "paused_step_budget"; break
            if status == "paused_step_budget":
                break
            before = np.asarray(metadata["cycle_start_state"])
            inv0,e0 = network.inventory_energy(before)
            inv1,e1 = network.inventory_energy(state)
            delta = state[network.process_size:] - before[network.process_size:]
            ns = network.ns
            feed,product,waste = delta[:ns],delta[ns:2*ns],delta[2*ns:3*ns]
            mass_residual = inv0+feed-product-waste-inv1
            energy_residual = e0+delta[3*ns]-delta[3*ns+1]-delta[3*ns+2]+delta[3*ns+3]-e1
            mass_relative = float(np.max(abs(mass_residual))/max(float(feed.sum()+product.sum()+waste.sum()),1e-6))
            energy_relative = abs(float(energy_residual))/max(float(np.sum(abs(delta[3*ns:3*ns+4]))),1.)
            closure = scaled_state_change(network,before,state)
            reporting_name = network.reporting_species
            reporting_index = network.reporting_index
            purity = float(product[reporting_index]/product.sum()) if product.sum()>1e-14 else None
            recovery = float(product[reporting_index]/feed[reporting_index]) if feed[reporting_index]>1e-14 else None
            previous = metadata["cycle_rows"][-1] if metadata["cycle_rows"] else None
            if previous:
                purity_change=0. if purity is None and previous.get("reporting_purity") is None else (abs(purity-previous["reporting_purity"]) if purity is not None and previous.get("reporting_purity") is not None else 1.)
                recovery_change=0. if recovery is None and previous.get("reporting_recovery") is None else (abs(recovery-previous["reporting_recovery"]) if recovery is not None and previous.get("reporting_recovery") is not None else 1.)
                output_change=max(purity_change,recovery_change)
            else:
                output_change=1.
            passed = (closure["state_scaled_max"] <= solver["css_tolerance"]
                      and closure["pressure_max_abs_bar"] <= solver["css_tolerance"]*case["high_pressure_bar"]
                      and closure["temperature_max_abs_k"] <= solver["css_tolerance"]*case["bed"]["ambient_k"]
                      and output_change <= solver["css_tolerance"]
                      and mass_relative <= solver["mass_balance_tolerance"]
                      and energy_relative <= solver["energy_balance_tolerance"]
                      and (recovery is None or 0<=recovery<=1+solver["mass_balance_tolerance"]))
            metadata["css_streak"] = metadata["css_streak"]+1 if passed else 0
            row = {"cycle":cycle+1,"h2_purity":purity if reporting_name == "H2" else None,
                   "h2_recovery":recovery if reporting_name == "H2" else None,
                   "reporting_species":reporting_name,
                   "reporting_purity":purity,"reporting_recovery":recovery,
                   "mean_external_feed_mol_s":float(feed.sum()/network.schedule["cycle_time_s"]),
                   "h2_product_mol_s":float(product[reporting_index]/network.schedule["cycle_time_s"]) if reporting_name == "H2" else None,
                   "reporting_product_mol_s":float(product[reporting_index]/network.schedule["cycle_time_s"]),
                   "feed_by_component_mol":feed.tolist(),"product_by_component_mol":product.tolist(),"tail_by_component_mol":waste.tolist(),
                   # The last four ledger entries are accumulated sensible-
                   # enthalpy/work terms.  Preserve them explicitly instead
                   # of forcing downstream bridges to reconstruct energy from
                   # component totals and an assumed temperature.  They are
                   # still model-ledger quantities (not a validated DWSIM
                   # property-package enthalpy) and therefore retain the
                   # existing fail-closed release policy.
                   "feed_enthalpy_j":float(delta[3*ns]),
                   "product_enthalpy_j":float(delta[3*ns+1]),
                   "tail_enthalpy_j":float(delta[3*ns+2]),
                   "ambient_heat_j":float(delta[3*ns+3]),
                   # The vacuum-work accumulator is the final ledger entry;
                   # the preceding entry is ambient-wall heat.  Keep both
                   # explicit so a downstream bridge cannot confuse the two
                   # energy terms.
                   "vacuum_work_j":float(delta[3*ns+4]),
                   "mass_balance_abs_mol":float(np.max(abs(mass_residual))),"mass_balance_relative":mass_relative,
                   "energy_balance_abs_j":abs(float(energy_residual)),"energy_balance_relative":energy_relative,
                   "energy_ledger_basis":case["bed"].get("energy_mode", "unspecified") + (
                       ": sensible plus Peng-Robinson departure internal-energy/enthalpy ledger plus vacuum work"
                       if case["bed"].get("pr_departure_energy", False)
                       else ": sensible enthalpy ledger plus vacuum work"
                   ),
                   "css_output_change":output_change,
                   **closure,"css_streak":metadata["css_streak"]}
            metadata["cycle_rows"].append(row)
            metadata["cycle_index"] = cycle+1; metadata["slot_index"] = 0
            metadata["cycle_start_state"] = state.tolist()
            metadata["wall_seconds_total"] = initial_wall_seconds+time.perf_counter()-started
            save_checkpoint(output,metadata,state)
            print(f"{case['case_id']} cycle={cycle+1} purity={purity} recovery={recovery} CSS={closure['state_scaled_max']:.3e} mass={mass_relative:.2e} energy={energy_relative:.2e}",flush=True)
            cycles_done += 1
            if noncyclic:
                status = "experiment_completed"; break
            if cycle+1>=solver["minimum_cycles"] and metadata["css_streak"]>=solver["css_consecutive_cycles"]:
                status = "css_converged"; break
            if cycle_budget is not None and cycles_done>=cycle_budget:
                status = "paused_cycle_budget"; break
        final = metadata["cycle_rows"][-1] if metadata["cycle_rows"] else None
        trajectory_screen = trajectory_domain_screen(input_screen, metadata.get("trajectory_extrema", {}), status)
        trajectory_screen.update(hashes)
        atomic_json(output/"trajectory_domain_screen.json", trajectory_screen)
        valid = status in {"css_converged", "experiment_completed"} and final is not None and final.get("reporting_purity") is not None
        result = {"schema":"dynamic-psa-result/v1","case_id":case["case_id"],"status":status,**hashes,
                  "number_of_beds":network.nb,"number_of_cycle_steps":network.schedule["number_of_cycle_steps"],
                  "reporting_species":network.reporting_species,
                  "schedule_contract":network.schedule_contract,
                  "completed_cycles":metadata["cycle_index"],"completed_slot_in_cycle":metadata["slot_index"],
                  "wall_seconds_total":metadata["wall_seconds_total"],"final_cycle":final,
                  "training_eligible":False,"css_training_candidate":valid,"resume_count":metadata["resume_count"],
                  "dataset_release_gate":"closed_pending_independent_resolution_and_source_validation_audit",
                  "trajectory_domain_sampled_within_documented":trajectory_screen["accepted_nodes_all_within_documented_domains"],
                  "trajectory_domain_sha256":digest(output/"trajectory_domain_screen.json"),
                  "state_css_converged":status=="css_converged",
                  "source_code_zip_sha256":source_zip_hash,
                  "energy_mode":network.bed.energy_mode,"equation_of_state":network.bed.equation_of_state,
                  "product_receiver_dynamic":network.receiver is not None,
                  "mesh":{"cell_centers_m":network.bed.cell_centers_m.tolist(),"cell_widths_m":network.bed.dz.tolist(),
                          "axial_dispersion_m2_s":network.bed.axial_dispersion_m2_s,
                          "axial_dispersion_model":network.bed.axial_dispersion_model,
                          "dispersion_particle_radius_m":network.bed.dispersion_particle_radius_m,
                          "layer_names":[m["name"] for m in network.bed.layer_indices]},
                  "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform()},
                  "claim":case["claim"]}
        atomic_json(result_path,result)
        write_csv(output/"cycles.csv",metadata["cycle_rows"])
        write_csv(output/"steps.csv",metadata["step_rows"])
        write_csv(output/"last_cycle_trace.csv",metadata["trace_rows"])
        (output/"active_step.json").unlink(missing_ok=True)
        return result
    except Exception as error:
        atomic_json(output/"failure.json",{**hashes,"error_type":type(error).__name__,"error":str(error),
                    "structured_details":getattr(error,"report",None),
                    "checkpoint_available":(output/"checkpoint.json").exists(),"training_eligible":False})
        raise
    finally:
        lock_path.unlink(missing_ok=True)
