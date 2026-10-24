"""Case-parallel driver; the same command supports local subprocesses/PBS.

No MPI dependency: independent cases offer coarser, more useful parallelism
than distributing a small stiff ODE across nodes. Each worker is a new Python
process and writes only its own output directory. No shared mutable globals.
"""
from __future__ import annotations
import os
# Limit BLAS oversubscription before importing any NumPy/SciPy module.
for variable in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[variable]="1"

import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import platform
from pathlib import Path
import subprocess
import sys

from .runner import atomic_json, digest, code_hash


def prepare_entry(entry,root):
    """Resolve/validate/create one output serially before dispatcher threads.

    Keep containment checks intact. Existing Windows directory identities are
    resolved a second time after creation; any discrepancy includes the exact
    paths in its error instead of an uninformative generic rejection.
    """
    root=root.resolve()
    case_path=(root/entry["case"]).resolve()
    output=(root/entry["output"]).resolve()
    if not case_path.is_relative_to(root) or not output.is_relative_to(root):
        raise ValueError(f"Manifest paths must stay within campaign root: root={str(root)!r}, case={str(case_path)!r}, output={str(output)!r}")
    if output==root or case_path==root:
        raise ValueError("A case/output must not be the campaign root itself")
    if digest(case_path)!=entry["input_sha256"]:
        raise ValueError(f"Frozen input hash differs: {case_path}")
    output.mkdir(parents=True,exist_ok=True)
    observed=output.resolve()
    if observed!=output or not observed.is_relative_to(root):
        raise ValueError(f"Output identity changed after creation: root={str(root)!r}, before={str(output)!r}, after={str(observed)!r}")
    return case_path,output


def execute(entry, root, expected_code_hash=None, cycle_budget=None, step_budget=None, prepared=None):
    if expected_code_hash is not None and code_hash()!=expected_code_hash:
        raise ValueError("Code changed while dispatching campaign; no mixed-version worker was launched")
    case_path,output=prepare_entry(entry,root) if prepared is None else prepared
    if digest(case_path)!=entry["input_sha256"]:
        raise ValueError(f"Frozen input hash differs: {case_path}")
    command=[sys.executable,"-u","-m","psa_dynamic","run","--case",str(case_path),"--output",str(output)]
    if cycle_budget is not None:
        command.extend(["--cycle-budget",str(cycle_budget)])
    if step_budget is not None:
        command.extend(["--step-budget",str(step_budget)])
    with (output/"worker.log").open("a",encoding="utf-8") as log:
        completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=False)
    result={"case":entry["case"],"output":entry["output"],"exit_code":completed.returncode}
    atomic_json(output/"worker_status.json",result)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--index",type=int)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--cycle-budget",type=int,help="Pause each worker after this invocation's completed cycles; does not change its frozen case")
    parser.add_argument("--step-budget",type=int)
    parser.add_argument("--require-python", help="Fail closed unless the dispatcher Python version exactly matches this value, e.g. 3.11.17")
    args=parser.parse_args()
    if args.require_python is not None and platform.python_version() != args.require_python:
        raise ValueError(f"Python version mismatch: required {args.require_python}, observed {platform.python_version()}")
    manifest=args.manifest.resolve()
    content=json.loads(manifest.read_text(encoding="utf-8"))
    if content.get("code_sha256")!=code_hash():
        raise ValueError("Campaign code hash differs from the frozen code; rebuild and re-test the campaign")
    entries=content["cases"]
    if len({item["output"] for item in entries})!=len(entries):
        raise ValueError("Duplicate output directories in campaign manifest")
    if any(value is not None and value<1 for value in (args.cycle_budget,args.step_budget)):
        raise ValueError("Execution budgets must be positive")
    root=manifest.parent
    canonical_outputs=[(root/item["output"]).resolve() for item in entries]
    if len(set(canonical_outputs))!=len(canonical_outputs):
        raise ValueError("Canonical duplicate output directories in campaign manifest; no worker launched")
    if args.index is not None:
        if args.index<0 or args.index>=len(entries):
            raise ValueError("Array index is outside the frozen manifest")
        prepared=prepare_entry(entries[args.index],root)
        result=execute(entries[args.index],root,content["code_sha256"],args.cycle_budget,args.step_budget,prepared)
        print(json.dumps(result));raise SystemExit(result["exit_code"])
    if args.workers<1:
        raise ValueError("workers must be positive")
    # Path resolution/creation happens on this thread, not concurrently in
    # Windows path fallback logic. One bad path does not stop safe cases.
    results=[None]*len(entries)
    plans={}
    seen={}
    for index,item in enumerate(entries):
        try:
            prepared=prepare_entry(item,root)
            if prepared[1] in seen:
                raise ValueError(f"Canonical output alias duplicates case {seen[prepared[1]]}")
            seen[prepared[1]]=index
            plans[index]=prepared
        except Exception as error:
            results[index]={"case":item["case"],"exit_code":1,"error_type":type(error).__name__,"error":str(error),"stage":"dispatch_preparation"}
    # Reject *all* entries in an aliased pair, not just the later writer.
    duplicates=[i for i,r in enumerate(results) if r is not None and "Canonical output alias" in r.get("error","")]
    if duplicates:
        raise ValueError("Canonical duplicate output directories in campaign manifest; no worker launched")
    # Threads only dispatch subprocesses; all numerical work is isolated.
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(execute,entries[i],root,content["code_sha256"],args.cycle_budget,args.step_budget,plan):i for i,plan in plans.items()}
        for future in as_completed(futures):
            index=futures[future]
            try:
                results[index]=future.result()
            except Exception as error:
                results[index]={"case":entries[index]["case"],"exit_code":1,"error_type":type(error).__name__,"error":str(error),"stage":"worker_dispatch"}
            atomic_json(root/"batch_progress.json",{"manifest_sha256":digest(manifest),"cases":results})
    atomic_json(root/"batch_status.json",{"manifest_sha256":digest(manifest),"cases":results})
    raise SystemExit(int(any(item["exit_code"]!=0 for item in results)))


if __name__=="__main__":
    main()
