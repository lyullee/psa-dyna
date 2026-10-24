import json
from pathlib import Path
import sys
import threading
import pytest
from psa_dynamic import batch
from psa_dynamic.runner import atomic_json,code_hash,digest


def test_containment_rejects_escape_and_root_without_external_writes(tmp_path):
    root=tmp_path/"campaign";root.mkdir()
    source=root/"case.json";atomic_json(source,{"case":"dummy"})
    entry={"case":"case.json","output":"../escape","input_sha256":digest(source)}
    with pytest.raises(ValueError,match="root=.*case=.*output="):
        batch.prepare_entry(entry,root)
    assert not (tmp_path/"escape").exists()
    entry["output"]="."
    with pytest.raises(ValueError,match="root itself"):
        batch.prepare_entry(entry,root)


def test_parallel_dispatch_resolves_paths_only_on_main_thread_and_isolates_bad_path(tmp_path,monkeypatch):
    source=tmp_path/"case.json";atomic_json(source,{"case":"dummy"})
    entries=[{"case":"case.json","output":out,"input_sha256":digest(source)} for out in ("good","../escape")]
    manifest=tmp_path/"manifest.json";atomic_json(manifest,{"cases":entries,"code_sha256":code_hash()})
    main=threading.get_ident()
    original=Path.resolve
    def serial_resolve(self,*args,**kwargs):
        assert threading.get_ident()==main,"Path resolution leaked into a dispatcher thread"
        return original(self,*args,**kwargs)
    monkeypatch.setattr(Path,"resolve",serial_resolve)
    def fake_execute(entry,root,expected,cycles,steps,prepared):
        assert prepared[0].is_absolute() and prepared[1].exists()
        return {"case":entry["case"],"exit_code":0}
    monkeypatch.setattr(batch,"execute",fake_execute)
    monkeypatch.setattr(sys,"argv",["batch","--manifest",str(manifest),"--workers","2"])
    with pytest.raises(SystemExit) as result:
        batch.main()
    assert result.value.code==1
    rows=json.loads((tmp_path/"batch_status.json").read_text())["cases"]
    assert [r["exit_code"] for r in rows]==[0,1]
    assert rows[1]["stage"]=="dispatch_preparation"


def test_canonical_output_alias_rejected_before_any_worker(tmp_path,monkeypatch):
    source=tmp_path/"case.json";atomic_json(source,{})
    entries=[{"case":"case.json","output":out,"input_sha256":digest(source)} for out in ("same","./same")]
    manifest=tmp_path/"manifest.json";atomic_json(manifest,{"cases":entries,"code_sha256":code_hash()})
    monkeypatch.setattr(sys,"argv",["batch","--manifest",str(manifest),"--workers","2"])
    monkeypatch.setattr(batch,"execute",lambda *a,**k:pytest.fail("No aliased writer may launch"))
    with pytest.raises(ValueError,match="Canonical duplicate"):
        batch.main()


def test_required_python_version_fails_closed_before_dispatch(tmp_path, monkeypatch):
    source = tmp_path / "case.json"
    atomic_json(source, {})
    manifest = tmp_path / "manifest.json"
    atomic_json(manifest, {"cases": [], "code_sha256": code_hash()})
    monkeypatch.setattr(sys, "argv", ["batch", "--manifest", str(manifest), "--require-python", "0.0.0"])
    with pytest.raises(ValueError, match="Python version mismatch"):
        batch.main()
