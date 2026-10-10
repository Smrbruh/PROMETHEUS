"""Reoptimize PROMETHEUS-01-v0 on h=1 mm, then freeze and audit at h=1/0.5.

This preregistered additional discretization study keeps r=6 mm and every other
canonical SIMP setting, including the 300 evaluated-design ceiling. It is not
the full h=2/1/0.5 fixed-field acceptance sequence and never grants that status.
Run sequentially after the canonical baseline to avoid competing CPU workloads.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import platform
import shutil
import sys
from time import perf_counter, process_time
import traceback

# Direct invocation from scripts/ must work without an editable installation.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import scipy
from prometheus.evaluation import Candidate, evaluate_candidate, verification_summary
from prometheus.geometry import make_l_mesh, bracket_pads, bracket_fixed_dofs, bracket_load_basis
from prometheus.models import Material, Benchmark
from prometheus.optimization import DensityFilter, solve_simp
from prometheus.run import jsonable, save_json, resource_checkpoint, peak_memory_bytes


def two_level_diagnostics(records):
    """Informative changes only; two levels cannot establish an error trend."""
    by_h = {row.get('mesh_h_mm'): row for row in records}
    if any(h not in by_h for h in (1., .5)):
        return {'status': 'resolution_unverified', 'reason': 'Two-level audit incomplete',
                'complete_v3_acceptance': False, 'diagnostic_thresholds_met': False}
    coarse, fine = by_h[1.], by_h[.5]
    if any(not row.get('solver_passed') or not row.get('independent', {}).get('solver_passed')
           for row in (coarse, fine)):
        return {'status': 'resolution_unverified', 'reason': 'A numerical solve check failed',
                'complete_v3_acceptance': False, 'diagnostic_thresholds_met': False}

    def change(first, second, name):
        return abs(first[name]-second[name])/max(abs(second[name]), 1e-12)

    values = {
        'q4_compliance_h1_to_h0_5': change(coarse, fine, 'compliance_n_mm'),
        'q4_pad_h1_to_h0_5': change(coarse, fine, 'max_pad_displacement_mm'),
        'p2_compliance_h1_to_h0_5': change(coarse['independent'], fine['independent'],
                                         'compliance_n_mm'),
        'p2_pad_h1_to_h0_5': change(coarse['independent'], fine['independent'],
                                  'max_pad_displacement_mm'),
        'finest_compliance_agreement': fine['independent']['q4_relative_compliance_difference'],
        'finest_pad_agreement': fine['independent']['q4_relative_pad_difference'],
    }
    criteria = (values['q4_compliance_h1_to_h0_5'] < .01,
                values['q4_pad_h1_to_h0_5'] < .02,
                values['p2_compliance_h1_to_h0_5'] < .01,
                values['p2_pad_h1_to_h0_5'] < .02,
                values['finest_compliance_agreement'] < .01,
                values['finest_pad_agreement'] < .02)
    return {**values, 'status': 'resolution_unverified',
            'diagnostic_thresholds_met': bool(all(criteria)),
            'complete_v3_acceptance': False,
            'reason': 'Only two frozen-field mesh levels; no third-level refinement trend',
            'physical_validation': 'not_assessed'}


def run_study(output):
    output = Path(output).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError('Output must be new or empty; preserve existing evidence')
    output.mkdir(parents=True, exist_ok=True)
    started, cpu = perf_counter(), process_time()
    material, benchmark = Material(), Benchmark()
    source_mesh = make_l_mesh(1.)
    pads = bracket_pads(source_mesh)
    initial = DensityFilter(source_mesh, pads, benchmark.filter_radius).initialize(
        benchmark.volume_fraction, benchmark.xmin)

    sources = sorted((ROOT/'prometheus').glob('*.py')) + [Path(__file__).resolve()]
    source_hashes = {}
    for path in sources:
        relative = path.relative_to(ROOT)
        source_hashes[relative.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        target = output/'source_snapshot'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    specifications = {}
    for name in ('MATHEMATICAL_FOUNDATION.md', 'BASELINE_EXPERIMENT.md',
                 'VALIDATION_PROTOCOL.md', 'EXPERIMENT_PROTOCOL.md'):
        path = ROOT/'research'/name
        specifications[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        target = output/'specification_snapshot'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    manifest = {
        'schema_version': 1, 'benchmark_id': 'PROMETHEUS-01-v0',
        'study_id': 'h1_reoptimization_two_level_frozen_audit',
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'command': sys.argv, 'material': asdict(material), 'benchmark': asdict(benchmark),
        'optimization_mesh_h_mm': 1., 'analysis_meshes_mm': [1., .5],
        'frozen_field_rule': 'Each child inherits its h1 parent physical density; no refilter/reoptimization during audit',
        'canonical_comparison_rule': 'Compare with h2-optimized B3 at common h0.5 and the same element/representation',
        'acceptance_scope': 'Two-level diagnostic only; full h2/1/0.5 V3 acceptance is not granted',
        'initialization': 'Constant free raw density from physical-volume bisection, preserved pads raw=physical=1',
        'units': {'length': 'mm', 'force': 'N', 'modulus': 'MPa', 'mass': 'kg', 'compliance': 'N mm'},
        'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
        'platform': platform.platform(), 'source_sha256': source_hashes,
        'specification_sha256': specifications,
        'threads': {key: os.environ.get(key) for key in
                    ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
        'numerical_backends': {'numpy': np.__config__.show(mode='dicts'),
                               'scipy': scipy.__config__.show(mode='dicts')},
        'resource_policy': '8 GiB measured peak checkpoint and P2 assembly guard; no OS hard limit',
        'parallelism_policy': 'Run after canonical baseline; no concurrent benchmark processes',
    }
    save_json(output/'manifest.json', manifest)
    np.savez_compressed(output/'source_mesh.npz', nodes=source_mesh.nodes,
                        cells=source_mesh.cells, cell_indices=source_mesh.cell_indices,
                        h_mm=np.array(source_mesh.h), pads=pads,
                        fixed_dofs=bracket_fixed_dofs(source_mesh),
                        unit_load_basis=bracket_load_basis(source_mesh))
    optimization_dir = output/'optimization_1'
    optimization_dir.mkdir()
    np.save(optimization_dir/'initial_raw.npy', initial)
    import json

    def callback(row):
        row['elapsed_wall_seconds'] = perf_counter()-optimization_started
        with (optimization_dir/'progress.jsonl').open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(jsonable(row), allow_nan=False)+'\n')
        if row['evaluation'] == 1 or row['evaluation'] % 25 == 0:
            print(f"h1 SIMP evaluation {row['evaluation']}: C={row['compliance']:.9g}, "
                  f"V={row['volume_fraction']:.9f}", flush=True)
        resource_checkpoint()

    optimization_started = perf_counter()
    try:
        optimum = solve_simp(source_mesh, material, benchmark, initial=initial, callback=callback)
    except Exception as exc:
        failure = {'status': 'failed', 'exception': type(exc).__name__, 'message': str(exc),
                   'traceback': traceback.format_exc(),
                   'history': getattr(exc, 'optimization_history', None),
                   'wall_seconds': perf_counter()-started,
                   'cpu_seconds': process_time()-cpu, 'peak_memory_bytes': peak_memory_bytes()}
        save_json(optimization_dir/'failure.json', failure)
        for attribute in ('last_evaluated_raw', 'failed_raw'):
            values = getattr(exc, attribute, None)
            if values is not None:
                np.save(optimization_dir/f'{attribute}.npy', values)
        save_json(output/'failure.json', failure)
        raise
    optimization_wall_seconds = perf_counter()-optimization_started
    for name in ('raw', 'physical', 'best_raw', 'best_physical'):
        np.save(optimization_dir/f'{name}.npy', getattr(optimum, name))
    save_json(optimization_dir/'history.json', optimum.history)
    keys = ['evaluation', 'compliance', 'volume', 'volume_fraction', 'volume_violation',
            'max_change', 'relative_compliance_spread', 'solve_seconds', 'elapsed_wall_seconds']
    with (optimization_dir/'history.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows({key: row[key] for key in keys} for row in optimum.history)
    optimization_summary = {
        'status': optimum.status, 'evaluations': optimum.evaluations,
        'wall_seconds': optimization_wall_seconds, 'solve_seconds': optimum.solve_seconds,
        'best_evaluation': optimum.best_evaluation,
        'final_compliance_n_mm': optimum.history[-1]['compliance'],
        'final_max_change': optimum.history[-1]['max_change'],
        'final_relative_compliance_spread': optimum.history[-1]['relative_compliance_spread'],
        'best_compliance_n_mm': optimum.history[optimum.best_evaluation-1]['compliance'],
        'raw_sha256': hashlib.sha256(optimum.raw.astype('<f8').tobytes()).hexdigest(),
        'physical_sha256': hashlib.sha256(optimum.physical.astype('<f8').tobytes()).hexdigest(),
    }
    save_json(optimization_dir/'summary.json', optimization_summary)
    candidate = Candidate('B3_h1_reoptimized', optimum.physical, 'density')
    candidate_dir = output/'candidates'/candidate.name
    candidate_dir.mkdir(parents=True)
    np.save(candidate_dir/'physical_source.npy', candidate.physical)
    records, failures = [], []
    for spacing in (1., .5):
        print(f'Frozen h1 design: Q4 and independent P2 audit at h={spacing:g} mm', flush=True)
        try:
            record, arrays = evaluate_candidate(candidate, source_mesh, material, benchmark,
                                                h=spacing, independent=True)
            record['process_peak_memory_bytes'] = resource_checkpoint()
            record['independent_verified'] = False
            record['acceptance_scope'] = 'two_level_diagnostic_only'
            np.savez_compressed(candidate_dir/f'h{spacing:g}_fields.npz', **arrays)
            save_json(candidate_dir/f'h{spacing:g}.json', record)
            records.append(record)
        except Exception as exc:
            failure = dict(getattr(exc, 'partial_evaluation', {}))
            failure.update({'status': 'failed', 'candidate': candidate.name, 'mesh_h_mm': spacing,
                            'exception': type(exc).__name__, 'message': str(exc),
                            'traceback': traceback.format_exc()})
            save_json(candidate_dir/f'h{spacing:g}_failure.json', failure)
            partial_arrays = getattr(exc, 'partial_arrays', None)
            if partial_arrays:
                np.savez_compressed(candidate_dir/f'h{spacing:g}_partial_fields.npz', **partial_arrays)
            failures.append(failure)
            records.append(failure)
            break
        finally:
            save_json(output/'evaluations_partial.json', records)
    # Common acceptance helper deliberately refuses this incomplete V3 sequence.
    acceptance = verification_summary(records)
    diagnostics = two_level_diagnostics(records)
    cost = {
        'optimization_factorizations': optimum.evaluations,
        'optimization_rhs_solves': optimum.evaluations,
        'evaluation_factorizations': sum(row.get('factorizations', 0) for row in records),
        'evaluation_rhs_solves': sum(row.get('rhs_solves', 0) for row in records),
        'failed_audit_attempts': len(failures),
        'wall_seconds': perf_counter()-started, 'cpu_seconds': process_time()-cpu,
        'peak_memory_bytes': peak_memory_bytes(),
        'failed_attempt_cost_note': 'Failed attempts remain in wall/CPU totals; unknown completed solve stages are not invented',
    }
    cost['completed_factorizations'] = cost['optimization_factorizations']+cost['evaluation_factorizations']
    cost['completed_rhs_solves'] = cost['optimization_rhs_solves']+cost['evaluation_rhs_solves']
    summary = {
        'manifest': 'manifest.json', 'study_id': manifest['study_id'],
        'optimization': optimization_summary, 'evaluations': records,
        'verification': acceptance, 'two_level_diagnostics': diagnostics,
        'cost': cost, 'failures': failures,
        'scientific_status': 'numerical_model_only_no_physical_validation',
        'conclusion_scope': 'Optimization mesh dependence screening; no mesh-independent or fully verified acceptance claim',
    }
    save_json(candidate_dir/'evaluations.json', records)
    save_json(candidate_dir/'verification.json', acceptance)
    save_json(output/'summary.json', summary)
    print(json.dumps({'output': str(output), 'optimization': optimization_summary,
                      'verification': acceptance, 'two_level_diagnostics': diagnostics,
                      'cost': cost}, indent=2), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='New or empty evidence directory')
    args = parser.parse_args()
    run_study(args.output)


if __name__ == '__main__':
    main()
