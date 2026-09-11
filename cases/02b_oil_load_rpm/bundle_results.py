"""Bundle compact completed evidence; leave every raw VTK trial on the workstation."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parent


def main():
    """Record executable identity and package final configurations, fields and logs."""
    settings=json.loads((ROOT/"settings.json").read_text())
    exe=Path(settings["executable"])
    rows=[]
    for pattern in ("scan_n*.json","fine_n*.json","low_n*.json"):
        rows.extend(json.loads(p.read_text()) for p in ROOT.glob(pattern))
    source=exe.parent.parent
    provenance=dict(host=platform.node(),run_root=str(ROOT),executable=str(exe),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        source_commit=subprocess.check_output(["git","-C",str(source),"rev-parse","HEAD"],text=True).strip(),
        source_modifications=subprocess.check_output(["git","-C",str(source),"status","--short","src"],text=True),
        equilibria=len(rows),trial_evaluations=len(list(ROOT.glob("runs/*/*/metrics.json"))),
        settings=settings)
    (ROOT/"provenance.json").write_text(json.dumps(provenance,indent=2))
    files={p for p in ROOT.iterdir() if p.is_file() and p.suffix!=".gz"}
    for row in rows:
        run=ROOT/row["dir"]
        files.update(run/name for name in ("config.txt","run.log","metrics.json","fields.npz"))
        files.add(run/"results/diagnostics.csv")
    with tarfile.open(ROOT/"oil_scan_compact.tar.gz","w:gz") as archive:
        for path in sorted(files):
            archive.add(path,arcname=path.relative_to(ROOT))
    print(json.dumps(provenance,indent=2))


if __name__=="__main__":
    main()
