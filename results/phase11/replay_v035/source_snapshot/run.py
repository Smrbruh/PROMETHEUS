"""Reproduce the baseline and preserve every result, failure and audit cost.

Run `python -m prometheus.run --help`. No adaptive selector is implemented.
"""
import os
# Set before numerical imports; documented fixed single-thread baseline.
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import ctypes
import hashlib
import json
from pathlib import Path
import platform
import shutil
import sys
from time import perf_counter, process_time
import traceback
import numpy as np
import scipy

from .models import Material, Benchmark
from .geometry import make_l_mesh, bracket_pads, bracket_fixed_dofs, bracket_load_basis
from .optimization import DensityFilter, solve_simp, simple_l_design, binary_by_rank, sensitivity_audit
from .evaluation import Candidate, evaluate_candidate, verification_summary


def jsonable(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {str(k):jsonable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [jsonable(v) for v in value]
    return value


def save_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(jsonable(value),indent=2,allow_nan=False)+'\n',encoding='utf-8')


def peak_memory_bytes():
    if os.name=='nt':
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[
                (name,ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize',
                'QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage',
                'QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
        counters=Counters(); counters.cb=ctypes.sizeof(counters)
        kernel=ctypes.windll.kernel32
        kernel.GetCurrentProcess.argtypes=[]
        kernel.GetCurrentProcess.restype=wintypes.HANDLE
        read_memory=ctypes.windll.psapi.GetProcessMemoryInfo
        read_memory.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
        read_memory.restype=wintypes.BOOL
        ok=read_memory(kernel.GetCurrentProcess(),ctypes.byref(counters),counters.cb)
        return int(counters.PeakWorkingSetSize) if ok else None
    try:
        import resource
        value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value if sys.platform=='darwin' else value*1024)
    except ImportError:return None


def resource_checkpoint():
    peak=peak_memory_bytes()
    if peak is not None and peak>8*1024**3:
        raise MemoryError('Measured peak working set exceeded 8 GiB checkpoint cap')
    return peak


def scientific_history(history):
    return [{k:v for k,v in row.items() if 'seconds' not in k} for row in history]


def run_benchmark(output,fraction=.4,repeats=2,meshes=(2.,1.,.5),audit_derivatives=False,
                  max_evaluations=300):
    output=Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Output must be new or empty; preserve existing evidence')
    output.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    material,benchmark=Material(),Benchmark(volume_fraction=fraction,max_evaluations=max_evaluations)
    mesh=make_l_mesh(2.)
    started,cpu=perf_counter(),process_time()
    files=sorted((root/'prometheus').glob('*.py'))
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    for p in files:
        target=output/'source_snapshot'/p.name
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    for name in ('MATHEMATICAL_FOUNDATION.md','BASELINE_EXPERIMENT.md','VALIDATION_PROTOCOL.md','EXPERIMENT_PROTOCOL.md'):
        target=output/'specification_snapshot'/name
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/'research'/name,target)
    np.savez_compressed(output/'source_mesh.npz',nodes=mesh.nodes,cells=mesh.cells,
                        cell_indices=mesh.cell_indices,h_mm=np.array(mesh.h),
                        pads=bracket_pads(mesh),fixed_dofs=bracket_fixed_dofs(mesh),
                        unit_load_basis=bracket_load_basis(mesh))
    initial_raw=DensityFilter(mesh,bracket_pads(mesh),benchmark.filter_radius).initialize(fraction,benchmark.xmin)
    np.save(output/'initial_raw.npy',initial_raw)
    manifest={'schema_version':1,'benchmark_id':'PROMETHEUS-01-v0',
              'research_date':'2026-10-10','started_utc':datetime.now(timezone.utc).isoformat(),
              'command':sys.argv,'material':asdict(material),'benchmark':asdict(benchmark),
              'meshes_mm':list(meshes),'repeats':repeats,
              'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
              'platform':platform.platform(),'source_sha256':hashes,
              'units':{'length':'mm','force':'N','modulus':'MPa','mass':'kg','compliance':'N mm'},
              'threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
              'numerical_backends':{'numpy':np.__config__.show(mode='dicts'),
                                    'scipy':scipy.__config__.show(mode='dicts')},
              'specification_sha256':{name:hashlib.sha256((root/'research'/name).read_bytes()).hexdigest()
                  for name in ('MATHEMATICAL_FOUNDATION.md','BASELINE_EXPERIMENT.md','VALIDATION_PROTOCOL.md','EXPERIMENT_PROTOCOL.md')},
              'resource_policy':'8 GiB peak checkpoint, conservative P2 assembly guard; not an OS hard limit'}
    save_json(output/'manifest.json',manifest)
    all_runs=[]
    for run in range(repeats):
        print(f'SIMP run {run+1}/{repeats}, volume={fraction:.2f}',flush=True)
        sub=output/f'optimization_{run+1}'
        sub.mkdir()
        progress=sub/'progress.jsonl'
        def callback(row):
            row['elapsed_wall_seconds']=perf_counter()-begin
            with progress.open('a',encoding='utf-8') as f:
                f.write(json.dumps(jsonable(row),allow_nan=False)+'\n')
            if row['evaluation']==1 or row['evaluation']%25==0:
                print(f"  evaluation {row['evaluation']}: C={row['compliance']:.8g}, V={row['volume_fraction']:.8f}",flush=True)
            resource_checkpoint()
        begin=perf_counter()
        try: opt=solve_simp(mesh,material,benchmark,initial=initial_raw,callback=callback)
        except Exception as exc:
            save_json(sub/'failure.json',{'exception':type(exc).__name__,'message':str(exc),
                'traceback':traceback.format_exc(),'history':getattr(exc,'optimization_history',None)})
            raise
        np.save(sub/'raw.npy',opt.raw);np.save(sub/'physical.npy',opt.physical)
        np.save(sub/'best_raw.npy',opt.best_raw);np.save(sub/'best_physical.npy',opt.best_physical)
        save_json(sub/'history.json',opt.history)
        keys=['evaluation','compliance','volume','volume_fraction','volume_violation',
              'max_change','relative_compliance_spread','solve_seconds']
        keys+=['elapsed_wall_seconds']
        with (sub/'history.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
            w.writerows({k:row[k] for k in keys} for row in opt.history)
        result={'status':opt.status,'evaluations':opt.evaluations,'wall_seconds':perf_counter()-begin,
                'solve_seconds':opt.solve_seconds,'best_evaluation':opt.best_evaluation,
                'final_compliance_n_mm':opt.history[-1]['compliance'],
                'final_max_change':opt.history[-1]['max_change'],
                'final_relative_compliance_spread':opt.history[-1]['relative_compliance_spread'],
                'best_compliance_n_mm':opt.history[opt.best_evaluation-1]['compliance'],
                'physical_sha256':hashlib.sha256(opt.physical.astype('<f8').tobytes()).hexdigest(),
                'raw_sha256':hashlib.sha256(opt.raw.astype('<f8').tobytes()).hexdigest()}
        save_json(sub/'summary.json',result); all_runs.append((opt,result))
    opt=all_runs[0][0]
    reproducibility={'repeated':repeats>1,'identical_arrays':None,'identical_numerical_history':None}
    if repeats>1:
        reproducibility.update(identical_arrays=all(np.array_equal(opt.raw,o.raw) and np.array_equal(opt.physical,o.physical) for o,_ in all_runs[1:]),
            identical_numerical_history=all(scientific_history(opt.history)==scientific_history(o.history) for o,_ in all_runs[1:]))
    filt=DensityFilter(mesh,bracket_pads(mesh),benchmark.filter_radius)
    initial=filt.physical(filt.initialize(fraction,benchmark.xmin))
    threshold=(opt.physical>=.5);threshold[bracket_pads(mesh)]=True
    candidates=[Candidate('B0',np.ones(len(mesh.cells)),'binary'),
                Candidate('B1',simple_l_design(mesh,fraction),'binary'),
                Candidate('B2',initial,'density'),Candidate('B3',opt.physical,'density'),
                Candidate('B3_binary',binary_by_rank(mesh,opt.physical,fraction),'binary'),
                Candidate('B3_threshold',threshold,'binary')]
    evaluations={};verifications={}; failures=[]
    for candidate in candidates:
        sub=output/'candidates'/candidate.name;sub.mkdir(parents=True)
        np.save(sub/'physical_source.npy',candidate.physical)
        records=[]
        for h in ((2.,) if candidate.name=='B2' else meshes):
            print(f'Audit {candidate.name}, h={h:g} mm',flush=True)
            try:
                record,arrays=evaluate_candidate(candidate,mesh,material,benchmark,h=h,
                                                independent=candidate.name!='B2')
                record['process_peak_memory_bytes']=resource_checkpoint()
                save_json(sub/f'h{h:g}.json',record)
                if arrays:np.savez_compressed(sub/f'h{h:g}_fields.npz',**arrays)
                records.append(record)
                if record['status']=='geometry_failed':break
            except Exception as exc:
                error=dict(getattr(exc,'partial_evaluation',{}))
                error.update({'candidate':candidate.name,'mesh_h_mm':h,'status':'failed',
                       'exception':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
                )
                save_json(sub/f'h{h:g}_failure.json',error);failures.append(error);records.append(error)
                partial_arrays=getattr(exc,'partial_arrays',None)
                if partial_arrays:np.savez_compressed(sub/f'h{h:g}_partial_fields.npz',**partial_arrays)
                break
        evaluations[candidate.name]=records
        verifications[candidate.name]=verification_summary(records)
        for r in records:
            r['independent_verified']=verifications[candidate.name]['passed']
            if r.get('status')!='failed':save_json(sub/f"h{r['mesh_h_mm']:g}.json",r)
        save_json(sub/'verification.json',verifications[candidate.name])
        save_json(sub/'evaluations.json',records)
        save_json(output/'evaluations_partial.json',evaluations)
    extra={}
    # Frozen nominal B3 ersatz dependence: no reoptimization or density change.
    from dataclasses import replace
    for ratio in (1e-8,1e-4):
        r,_=evaluate_candidate(candidates[3],mesh,material,replace(benchmark,emin_ratio=ratio),
                                h=2.,independent=False,include_stress=False)
        extra[str(ratio)]=r
    derivatives={}
    if audit_derivatives:
        for signed in (False,True):
            print(f'Sensitivity finite differences, signed={signed}',flush=True)
            audit=sensitivity_audit(mesh,material,benchmark,signed=signed)
            name='signed' if signed else 'well_scaled'
            np.savez_compressed(output/f'sensitivity_{name}.npz',**audit)
            derivatives[name]=jsonable(audit)
        save_json(output/'sensitivity.json',derivatives)
    totals={'optimization_factorizations':sum(o.evaluations for o,_ in all_runs),
            'optimization_rhs_solves':sum(o.evaluations for o,_ in all_runs),
            'evaluation_factorizations':sum(r.get('factorizations',0) for rr in evaluations.values() for r in rr),
            'evaluation_rhs_solves':sum(r.get('rhs_solves',0) for rr in evaluations.values() for r in rr),
            'ersatz_factorizations':2,'ersatz_rhs_solves':4,
            'sensitivity_factorizations':202 if audit_derivatives else 0,
            'sensitivity_rhs_solves':202 if audit_derivatives else 0,
            'failed_audit_attempts':len(failures),
            'wall_seconds':perf_counter()-started,'cpu_seconds':process_time()-cpu,
            'peak_memory_bytes':peak_memory_bytes()}
    totals['completed_factorizations']=sum(totals[k] for k in ('optimization_factorizations','evaluation_factorizations','ersatz_factorizations','sensitivity_factorizations'))
    totals['completed_rhs_solves']=sum(totals[k] for k in ('optimization_rhs_solves','evaluation_rhs_solves','ersatz_rhs_solves','sensitivity_rhs_solves'))
    totals['failed_attempt_cost_note']='Failed attempts remain charged in wall/CPU time; factorization stage counts are unknown for failures.'
    summary={'manifest':'manifest.json','optimization':[r for _,r in all_runs],
             'reproducibility':reproducibility,'evaluations':evaluations,'verification':verifications,
             'ersatz_sensitivity':extra,'cost':totals,'failures':failures,
             'scientific_status':'numerical_model_only_no_physical_validation'}
    save_json(output/'summary.json',summary)
    print(json.dumps({'output':str(output),'optimization':opt.status,'evaluations':opt.evaluations,
                      'reproducibility':reproducibility,'verification':{k:v['status'] for k,v in verifications.items()},
                      'cost':totals},indent=2),flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--volume',type=float,default=.4)
    parser.add_argument('--repeats',type=int,choices=(1,2),default=2)
    parser.add_argument('--meshes',nargs='+',type=float,default=[2.,1.,.5])
    parser.add_argument('--audit-derivatives',action='store_true')
    parser.add_argument('--max-evaluations',type=int,default=300,help='Overrides label a pilot, not canonical completed baseline')
    args=parser.parse_args()
    run_benchmark(args.output,args.volume,args.repeats,tuple(args.meshes),args.audit_derivatives,args.max_evaluations)


if __name__=='__main__':main()
