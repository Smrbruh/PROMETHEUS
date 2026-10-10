"""Replay documented Phase 1 jobs into new directories and compare science."""
import os
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import json
import subprocess
import sys
from time import perf_counter
import hashlib
import argparse
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'phase11'


def strip_times(value):
    if isinstance(value,dict):return {k:strip_times(v) for k,v in value.items()
        if 'seconds' not in k and 'memory' not in k}
    if isinstance(value,list):return [strip_times(v) for v in value]
    return value


def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default=str(OUT))
    args=parser.parse_args();OUT=Path(args.output).resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    jobs=[('phase1_v040_verified','replay_v040',['-m','prometheus.run','--volume','.40','--repeats','2','--audit-derivatives']),
          ('phase1_v035','replay_v035',['-m','prometheus.run','--volume','.35','--repeats','1']),
          ('phase1_v045','replay_v045',['-m','prometheus.run','--volume','.45','--repeats','1']),
          ('phase1_h1_reoptimization','replay_h1',['scripts/mesh_dependence.py'])]
    comparisons=[]
    for old_name,new_name,args in jobs:
        old=ROOT/'results'/old_name;new=OUT/new_name
        if not (new/'summary.json').exists():
            if new.exists():raise FileExistsError(f'Incomplete output retained; choose a new run identity: {new}')
            print(f'Reproduction {new_name}',flush=True)
            with (OUT/f'{new_name}.log').open('x',encoding='utf-8') as log:
                result=subprocess.run([sys.executable,'-B',*args,'--output',str(new)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:raise RuntimeError(f'{new_name} failed; log retained')
        a=json.loads((old/'summary.json').read_text());b=json.loads((new/'summary.json').read_text())
        matched={key:strip_times(a[key])==strip_times(b[key]) for key in ['optimization','evaluations','verification','failures']}
        if 'reproducibility' in a:matched['reproducibility']=a['reproducibility']==b['reproducibility']
        arrays=[]
        for source in old.rglob('*.npy'):
            target=new/source.relative_to(old)
            if target.exists():arrays.append(bool(np.array_equal(np.load(source),np.load(target))))
        matched['all_npy_arrays_bit_identical']=bool(arrays) and all(arrays)
        comparisons.append({'original':old_name,'replay':new_name,'matches':matched,'arrays_compared':len(arrays),
            'passed':all(matched.values()),'replay_cost':b['cost'],
            'original_summary_sha256':hashlib.sha256((old/'summary.json').read_bytes()).hexdigest()})
        print(f'{new_name}: scientific match={comparisons[-1]["passed"]}',flush=True)
    validation=OUT/'validation_replay.json'
    if not validation.exists():
        subprocess.run([sys.executable,'-B','scripts/validation_evidence.py','--output',str(validation)],cwd=ROOT,check=True)
    report={'passed':all(c['passed'] for c in comparisons),'comparisons':comparisons,
            'validation_replay':'validation_replay.json','note':'Fresh process executions; exact numerical/history/label/array comparison excludes timing and process memory.'}
    target=OUT/'reproduction.json'
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if not report['passed']:raise AssertionError('Scientific reproduction mismatch; retained evidence')


if __name__=='__main__':main()
