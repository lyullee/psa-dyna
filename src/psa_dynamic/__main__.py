"""One entry point on Windows and Linux: python -m psa_dynamic ..."""
import argparse
import json
from pathlib import Path
from .cases import guderian_2021_n2_reference_case, s015_case, wang_case
from .runner import atomic_json, run_case
from .qualification import case_screen


def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="command",required=True)
    build=sub.add_parser("make-case")
    build.add_argument("--source",choices=["wang","s015","guderian"],required=True)
    build.add_argument("--schedule",required=True)
    build.add_argument("--cells",type=int,default=4)
    build.add_argument("--experiment-id",type=int,default=13,
                       help="Guderian workbook experiment ID (13, 14 or 15)")
    build.add_argument("--boundary-diagnostics", action=argparse.BooleanOptionalAction,
                       default=None,
                       help="Enable/disable detailed boundary quadrature diagnostics")
    build.add_argument("--output",type=Path,required=True)
    run=sub.add_parser("run")
    run.add_argument("--case",type=Path,required=True)
    run.add_argument("--output",type=Path,required=True)
    run.add_argument("--cycle-budget",type=int)
    run.add_argument("--step-budget",type=int)
    audit=sub.add_parser("audit-case",help="Report engineering evidence gaps without running a simulation")
    audit.add_argument("--case",type=Path,required=True)
    audit.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.command=="make-case":
        if args.source == "guderian":
            if args.schedule != "guderian_two_bed_6step_source":
                parser.error("Guderian cases require --schedule guderian_two_bed_6step_source")
            case = guderian_2021_n2_reference_case(args.experiment_id, args.cells)
        else:
            case=(wang_case if args.source=="wang" else s015_case)(args.schedule,args.cells)
        if args.boundary_diagnostics is not None:
            case["solver"]["record_boundary_diagnostics"] = args.boundary_diagnostics
        args.output.parent.mkdir(parents=True,exist_ok=True)
        atomic_json(args.output,case)
        print(args.output)
    elif args.command=="audit-case":
        report=case_screen(json.loads(args.case.read_text(encoding="utf-8-sig")))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        atomic_json(args.output,report)
        print(json.dumps(report,indent=2))
    else:
        result=run_case(args.case,args.output,cycle_budget=args.cycle_budget,step_budget=args.step_budget)
        print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
