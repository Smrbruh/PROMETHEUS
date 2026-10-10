"""Add-only, monitored fixed-geometry refinement; no Phase 1 writes."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import copy
import csv
import hashlib
import json
import subprocess
import time
import traceback
from dataclasses import asdict
import numpy as np
from prometheus.models import Material,Benchmark
from prometheus.geometry import make_l_mesh
from prometheus.evaluation import Candidate,evaluate_candidate
from prometheus.reference_labels import reference_label
from prometheus.audit_responses import project_response,resource_estimate
from prometheus.run import peak_memory_bytes

ANGLES=(-20.,0.,20.)
LEVELS=(2.,1.,.5,.25,.125,.0625)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_json(path,data):
    with Path(path).open('x',encoding='utf-8') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n')


def relative(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def skeleton(entry,h):
    source=make_l_mesh(2.)
    candidate=Candidate(entry['candidate_name'],np.load(ROOT/entry['array_path']),'binary')
    mat=Material();benchmark=Benchmark(volume_fraction=entry['volume_setting'])
    volume=mat.thickness*float(source.areas@candidate.physical)
    return {'candidate':candidate.name,'hash':candidate.content_hash,'representation':'binary',
        'source_h_mm':2.,'mesh_h_mm':h,'source_geometry_hash':hashlib.sha256(source.nodes.tobytes()+source.cells.tobytes()).hexdigest(),
        'problem_signature':hashlib.sha256(json.dumps({'material':asdict(mat),'benchmark':asdict(benchmark)},sort_keys=True).encode()).hexdigest(),
        'volume_fraction':volume/(mat.thickness*sum(source.areas)),
        'volume_budget_fraction':benchmark.volume_fraction,'volume_mm3':volume,'mass_kg':volume*mat.density,
        'volume_feasible':bool(volume<=benchmark.volume_fraction*32000+1e-6*32000),
        'pads_preserved':entry['geometry']['pads_preserved'],'status':'pending','solver_passed':False,
        'factorizations':0,'rhs_solves':0}


def worker(request_path):
    request=json.loads(Path(request_path).read_text(encoding='utf-8'))
    entry=request['entry'];h=request['mesh_h_mm'];out=Path(request_path).parent
    array_digest=digest(ROOT/entry['array_path'])
    expected=entry.get('provenance',{}).get('array_file_sha256',entry['source_input_sha256'])
    if array_digest!=expected:raise ValueError('Input digest changed')
    start=time.perf_counter();cpu=time.process_time();arrays={}
    try:
        record,arrays=evaluate_candidate(Candidate(entry['candidate_name'],np.load(ROOT/entry['array_path']),'binary'),
            make_l_mesh(2.),Material(),Benchmark(volume_fraction=entry['volume_setting']),h,True,False)
    except Exception as exc:
        record=getattr(exc,'partial_evaluation',skeleton(entry,h))
        arrays=getattr(exc,'partial_arrays',{})
        if record['status']=='pending':record['status']='solver_exception'
        record['exception']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    record.update(attempt_wall_seconds=time.perf_counter()-start,attempt_cpu_seconds=time.process_time()-cpu,
        process_peak_memory_bytes=peak_memory_bytes(),source_archive=entry['array_path'],source_input_sha256=entry['source_input_sha256'],
        analyzed_array_sha256=array_digest)
    # Preserve a completed numerical record before potentially expensive compression.
    write_json(out/'record.json',record)
    if arrays:
        with (out/'fields.npz').open('xb') as f:np.savez_compressed(f,**arrays)
    write_json(out/'worker_completion.json',{'wall_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu,
        'peak_memory_bytes':peak_memory_bytes(),'fields_saved':bool(arrays),'record_sha256':digest(out/'record.json')})


def working_set(process):
    if os.name!='nt':return None
    import ctypes
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[(k,ctypes.c_size_t) for k in
            ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage',
             'QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
    fn=ctypes.WinDLL('psapi',use_last_error=True).GetProcessMemoryInfo
    fn.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD];fn.restype=wintypes.BOOL
    counters=Counters();counters.cb=ctypes.sizeof(counters)
    return int(counters.WorkingSetSize) if fn(wintypes.HANDLE(int(process._handle)),ctypes.byref(counters),counters.cb) else None


def attempt(entry,h,out):
    out.mkdir(parents=True,exist_ok=False)
    estimate=resource_estimate(make_l_mesh(2.),np.load(ROOT/entry['array_path']).astype(bool),h)
    write_json(out/'request.json',{'entry':entry,'mesh_h_mm':h,'preflight':estimate,
        'protocol_sha256':digest(ROOT/'research/REFERENCE_LABEL_PROTOCOL.md'),
        'code_sha256':{relative(p):digest(p) for p in [Path(__file__),*sorted((ROOT/'prometheus').glob('*.py'))]},
        'policy':{'max_wall_seconds':600,'max_observed_working_set_bytes':8*1024**3,'max_p2_dofs':1500000}})
    start=time.perf_counter();peak=0;stop=None
    if not estimate['allowed']:
        record=skeleton(entry,h);record.update(status='resource_preflight_deferred',resource_preflight=estimate)
        write_json(out/'record.json',record)
        cost={'catalog_key':entry['catalog_key'],'mesh_h_mm':h,'wall_seconds':time.perf_counter()-start,
              'factorizations':0,'rhs_solves':0,'status':record['status'],'preflight':estimate,'peak_working_set_bytes':0}
        write_json(out/'cost.json',cost)
        return record,{},cost
    print(f"Refine {entry['catalog_key']} h={h:g}; estimated P2 DOFs {estimate['p2_dofs']:,}",flush=True)
    with (out/'worker.log').open('x',encoding='utf-8') as log:
        process=subprocess.Popen([sys.executable,'-B',str(Path(__file__)),'--worker',str(out/'request.json')],
            cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        while process.poll() is None:
            memory=working_set(process)
            if memory is not None:peak=max(peak,memory)
            if time.perf_counter()-start>600:stop='time_limit'
            elif memory is not None and memory>8*1024**3:stop='memory_limit'
            if stop:
                process.terminate();process.wait();break
            time.sleep(.25)
    record_path=out/'record.json';complete=(out/'worker_completion.json').exists()
    record=json.loads(record_path.read_text()) if record_path.exists() else skeleton(entry,h)
    if not complete:
        record=copy.deepcopy(record)
        record.update(status=stop or 'worker_incomplete',solver_passed=False,
            termination={'reason':stop or 'child_did_not_complete','exit_code':process.returncode},
            completed_solve_counts_known=record_path.exists())
        write_json(out/'terminal_record.json',record)
    completion=json.loads((out/'worker_completion.json').read_text()) if complete else {}
    arrays={}
    if complete and (out/'fields.npz').exists():
        with np.load(out/'fields.npz') as fields:arrays={k:fields[k] for k in fields.files}
    cost={'catalog_key':entry['catalog_key'],'mesh_h_mm':h,'wall_seconds':time.perf_counter()-start,
        'worker_cpu_seconds':completion.get('cpu_seconds'),'worker_wall_seconds':completion.get('wall_seconds'),
        'factorizations':record.get('factorizations') if complete or record_path.exists() else None,
        'rhs_solves':record.get('rhs_solves') if complete or record_path.exists() else None,
        'status':record['status'],'preflight':estimate,'peak_working_set_bytes':max(peak,completion.get('peak_memory_bytes') or 0),
        'record_path':relative(out/('terminal_record.json' if not complete else 'record.json')),
        'field_path':relative(out/'fields.npz') if complete and arrays else None,'postprocessing_complete':complete}
    write_json(out/'cost.json',cost)
    print(f"  {record['status']}; {cost['wall_seconds']:.2f}s; peak {cost['peak_working_set_bytes']/1024**3:.3f} GiB",flush=True)
    return record,arrays,cost


def main(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    catalog_path=ROOT/'results/phase11/geometry/catalog.json'
    catalog=json.loads(catalog_path.read_text())['entries'];start=time.perf_counter();cpu=time.process_time()
    write_json(out/'plan.json',{'catalogue_sha256':digest(catalog_path),'protocol_sha256':digest(ROOT/'research/REFERENCE_LABEL_PROTOCOL.md'),
        'levels_mm':LEVELS,'response_angles_degrees':ANGLES,
        'clarification':'Extend a connected binary if any declared response angle remains numerically unresolved; historical nominal acceptance does not shortcut added-angle checks. Stop only when all three angle gates pass, an attempt fails/resources defer, or planned levels end.',
        'created_before_refinements':True})
    labels=[];rows=[];costs=[];projection_seconds=0.;imports=[]
    for entry in catalog:
        histories={angle:[] for angle in ANGLES}
        if entry['history_path']:
            history_path=ROOT/entry['history_path'];old=json.loads(history_path.read_text())
            for record in old:
                fields_path=history_path.parent/f"h{record['mesh_h_mm']:g}_fields.npz";arrays={}
                if fields_path.exists():
                    with np.load(fields_path) as fields:arrays={k:fields[k] for k in fields.files}
                tick=time.perf_counter()
                for angle in ANGLES:histories[angle].append(project_response(record,arrays,angle))
                projection_seconds+=time.perf_counter()-tick
                imports.append({'catalog_key':entry['catalog_key'],'mesh_h_mm':record['mesh_h_mm'],
                    'record_archive':relative(history_path),'record_archive_sha256':digest(history_path),
                    'field_archive':relative(fields_path) if arrays else None,
                    'field_archive_sha256':digest(fields_path) if arrays else None,
                    'historical_wall_seconds':record.get('wall_seconds'),'new_factorizations':0,'new_rhs_solves':0})
        def current_labels():
            return {a:reference_label(histories[a],historical_label=entry.get('historical_label'),geometry_audit=entry['geometry']) for a in ANGLES}
        connected=entry['geometry']['geometry_valid']
        if connected:
            done={r['mesh_h_mm'] for r in histories[0.]}
            for h in LEVELS:
                if h in done:continue
                if len(histories[0.])>=3 and all(v['passed'] for v in current_labels().values()):break
                record,arrays,cost=attempt(entry,h,out/'attempts'/entry['catalog_key']/f'h{h:g}')
                costs.append(cost);tick=time.perf_counter()
                for angle in ANGLES:histories[angle].append(project_response(record,arrays,angle))
                projection_seconds+=time.perf_counter()-tick
                del arrays
                if record['status']!='evaluated':break
        final=current_labels()
        for angle in ANGLES:
            labels.append({'catalog_key':entry['catalog_key'],'candidate_name':entry['candidate_name'],
                'volume_setting':entry['volume_setting'],'angle_degrees':angle,'label':final[angle]})
            previous=None
            for record in histories[angle]:
                part=record.get('independent',{});C=record.get('compliance_n_mm');P=part.get('compliance_n_mm')
                def change(a,b):return abs(a-b)/abs(b) if a is not None and b not in (None,0) else None
                rows.append({'catalog_key':entry['catalog_key'],'candidate_name':entry['candidate_name'],
                    'volume_setting':entry['volume_setting'],'angle_degrees':angle,'mesh_h_mm':record['mesh_h_mm'],
                    'status':record['status'],'q4_compliance_n_mm':C,'p2_compliance_n_mm':P,
                    'q4_pad_mm':record.get('max_pad_displacement_mm'),'p2_pad_mm':part.get('max_pad_displacement_mm'),
                    'q4_compliance_change':change(previous.get('compliance_n_mm'),C) if previous else None,
                    'p2_compliance_change':change(previous.get('independent',{}).get('compliance_n_mm'),P) if previous else None,
                    'q4_pad_change':change(previous.get('max_pad_displacement_mm'),record.get('max_pad_displacement_mm')) if previous else None,
                    'p2_pad_change':change(previous.get('independent',{}).get('max_pad_displacement_mm'),part.get('max_pad_displacement_mm')) if previous else None,
                    'same_h_compliance_gap':change(C,P),'same_h_pad_gap':change(record.get('max_pad_displacement_mm'),part.get('max_pad_displacement_mm')),
                    'q4_solver_passed':record.get('solver_passed'), 'p2_solver_passed':part.get('solver_passed'),
                    'candidate_hash':record['hash'],'problem_signature':record['problem_signature'],
                    'volume_fraction':record['volume_fraction'],'components':entry['geometry']['edge_component_count'],
                    'min_axis_run_mm':entry['geometry']['minimum_axis_run_mm'],'local_contacts':entry['geometry']['local_point_contact_count']})
                previous=record
        candidate_out=out/'candidates'/entry['catalog_key'];candidate_out.mkdir(parents=True)
        write_json(candidate_out/'responses.json',{str(a):histories[a] for a in ANGLES})
        write_json(candidate_out/'labels.json',{str(a):final[a] for a in ANGLES})
        print(entry['catalog_key']+': '+', '.join(f"{a:g}deg={v['classification']}" for a,v in final.items()),flush=True)
    write_json(out/'labels.json',labels);write_json(out/'responses.json',rows)
    with (out/'responses.csv').open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write_json(out/'imported_evidence.json',imports)
    write_json(out/'costs.json',{'attempts':costs,'projection_wall_seconds':projection_seconds,
        'orchestrator_wall_seconds':time.perf_counter()-start,'parent_cpu_seconds':time.process_time()-cpu,
        'completed_factorizations':sum(c['factorizations'] or 0 for c in costs),
        'completed_rhs_solves':sum(c['rhs_solves'] or 0 for c in costs),
        'has_unknown_solve_counts':any(c['factorizations'] is None for c in costs)})
    with (out/'costs.csv').open('x',newline='',encoding='utf-8') as f:
        keys=['catalog_key','mesh_h_mm','status','wall_seconds','worker_cpu_seconds','factorizations','rhs_solves','peak_working_set_bytes']
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(costs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--worker');p.add_argument('--output',default='results/phase11/refinement')
    args=p.parse_args()
    if args.worker:worker(args.worker)
    else:main(args.output)
