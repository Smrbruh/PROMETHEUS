"""Capture/check Phase 1 file hashes; never writes to a Phase 1 artifact."""
from pathlib import Path
import argparse
import hashlib
import json
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'results'/'phase11'
OLD_CORE=['__init__','envelope','evaluation','fem','geometry','models','optimization','run','verification']
OLD_TESTS=['envelope','evaluation','fem','optimization','run','verification']
OLD_SCRIPTS=['check_phase1_evidence','mesh_dependence','report_evidence','validation_evidence']
NEW_DOCS={'REFERENCE_LABEL_PROTOCOL.md','NUMERICAL_TRUTH_AUDIT.md','BINARY_GEOMETRY_AUDIT.md','PHASE1_1_DECISION.md'}


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['capture','check'])
    p.add_argument('--output',default=None);a=p.parse_args()
    DEST.mkdir(parents=True,exist_ok=True)
    baseline=DEST/'preservation_before.json'
    if a.mode=='capture':
        if baseline.exists():raise FileExistsError(baseline)
        files=[ROOT/'prometheus'/f'{s}.py' for s in OLD_CORE]
        files += [ROOT/'tests'/f'test_{s}.py' for s in OLD_TESTS]
        files += [ROOT/'scripts'/f'{s}.py' for s in OLD_SCRIPTS]
        files += [ROOT/s for s in ['README.md','pyproject.toml','requirements-tested.txt','.gitignore']]
        for folder in ['research','results']:
            files += [f for f in (ROOT/folder).rglob('*') if f.is_file()
                      and '__pycache__' not in f.parts and DEST not in f.parents
                      and not (folder=='research' and f.name in NEW_DOCS)]
        evidence={'captured_utc':datetime.now(timezone.utc).isoformat(),'scope':'All preexisting source, tests, reports and archived evidence; runtime caches excluded',
                  'files':{f.relative_to(ROOT).as_posix():{'bytes':f.stat().st_size,'sha256':digest(f)} for f in sorted(set(files))}}
        baseline.write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
        print({'captured_files':len(evidence['files']),'path':str(baseline)})
    else:
        evidence=json.loads(baseline.read_text(encoding='utf-8'))
        failures=[r for r,v in evidence['files'].items() if not (ROOT/r).exists() or digest(ROOT/r)!=v['sha256']]
        report={'checked_utc':datetime.now(timezone.utc).isoformat(),'baseline_sha256':digest(baseline),'files_checked':len(evidence['files']),
                'unchanged':not failures,'changed_or_missing':failures}
        output=Path(a.output) if a.output else DEST/'preservation_after.json'
        if output.exists():raise FileExistsError(output)
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(report,indent=2))
        if failures:raise AssertionError('Original evidence changed')


if __name__=='__main__':main()
