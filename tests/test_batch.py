"""One malformed case must not stop other independent numerical workers."""
import json
from pathlib import Path
import subprocess
import sys

from psa_dynamic.cases import wang_case
from psa_dynamic.runner import atomic_json,code_hash,digest


def test_local_batch_isolates_case_failures(tmp_path):
    good=wang_case("closed",2);good["schedule"]["duration_s"]=.1
    bad=wang_case("closed",2);bad["low_pressure_bar"]=-1.
    entries=[]
    for name,case in [("good",good),("bad",bad)]:
        path=tmp_path/(name+".json");atomic_json(path,case)
        entries.append({"case":path.name,"output":name,"input_sha256":digest(path)})
    manifest=tmp_path/"manifest.json"
    atomic_json(manifest,{"cases":entries,"code_sha256":code_hash()})
    run=subprocess.run([sys.executable,"-m","psa_dynamic.batch","--manifest",str(manifest),"--workers","2"],capture_output=True,text=True)
    assert run.returncode==1
    status=json.loads((tmp_path/"batch_status.json").read_text())
    assert [r["exit_code"] for r in status["cases"]]==[0,1]
    assert json.loads((tmp_path/"good"/"result.json").read_text())["status"]=="experiment_completed"
    assert (tmp_path/"bad"/"worker.log").exists()
