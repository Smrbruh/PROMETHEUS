"""Append-only kinematics from successful repair Q4 fields after a P2 failure.

This does not rerun the catalogue, solve mechanics, or change an acceptance
label. The repaired geometry, Q4 result and failed P2 status remain distinct.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
from prometheus.kinematic_audit import audit_q4_kinematics
from scripts.phase11_kinematics import digest, relative, write_new, flatten_summary


def run(output, attempt):
    started, cpu = perf_counter(), process_time()
    output, attempt = Path(output).resolve(), Path(attempt).resolve()
    output.mkdir(parents=True, exist_ok=False)
    catalogue_path = ROOT/'results/phase11/geometry/catalog.json'
    source_record_path, fields_path = attempt/'record.json', attempt/'fields.npz'
    catalogue = json.loads(catalogue_path.read_text(encoding='utf-8'))
    entry = next(row for row in catalogue['entries'] if row['catalog_key'] == 'v035_repair')
    record = json.loads(source_record_path.read_text(encoding='utf-8'))
    if not record.get('solver_passed') or not record.get('diagnostics', {}).get('passed'):
        raise ValueError('Repair kinematics requires successful recorded Q4 algebraic checks')
    if record['candidate'] != entry['candidate_name'] or record['representation'] != 'binary':
        raise ValueError('Repair candidate identity does not match the catalogue')
    array_path = ROOT/entry['array_path']
    physical = np.load(array_path, allow_pickle=False)
    field_hash = hashlib.sha256(physical.astype('<f8').tobytes()+b'binary').hexdigest()
    if field_hash != record['hash'] or digest(array_path) != record['analyzed_array_sha256']:
        raise ValueError('Repair source array/content hash mismatch')
    input_paths = [catalogue_path, source_record_path, fields_path, array_path,
                   ROOT/entry['source_mesh_path']]
    input_hashes = {relative(path): digest(path) for path in input_paths}
    code_paths = [Path(__file__).resolve(), ROOT/'scripts/phase11_kinematics.py',
                  ROOT/'prometheus/kinematic_audit.py', ROOT/'prometheus/fem.py',
                  ROOT/'research/REFERENCE_LABEL_PROTOCOL.md']
    code_hashes = {relative(path): digest(path) for path in code_paths}
    with np.load(fields_path, allow_pickle=False) as archive:
        nodes, cells, basis, rho = (archive[name] for name in ('nodes', 'cells', 'u_basis', 'rho'))
    if basis.shape != (2*len(nodes), 2) or not np.all(rho == 1):
        raise ValueError('Expected two Q4 displacement bases on material-only binary cells')
    x, y = nodes.T
    pad = (((x >= -1e-10) & (x <= 40+1e-10) & (y >= 70-1e-10) & (y <= 80+1e-10))
           | ((x >= 110-1e-10) & (x <= 120+1e-10) & (y >= 10-1e-10) & (y <= 30+1e-10)))
    rows = []
    for angle in (-20., 0., 20.):
        theta = np.deg2rad(angle)
        u = (np.sin(theta)*basis[:, 0]+np.cos(theta)*basis[:, 1]).reshape(-1, 2)
        summary, _ = audit_q4_kinematics(nodes, cells, u)
        magnitude = np.linalg.norm(u, axis=1)
        maximum, pad_maximum = float(magnitude.max()), float(magnitude[pad].max())
        rows.append({'catalog_key': entry['catalog_key'], 'candidate_name': entry['candidate_name'],
            'catalog_identity_sha256': entry['catalog_identity_sha256'],
            'geometry_sha256': entry['geometry_sha256'], 'candidate_content_hash': record['hash'],
            'volume_setting': .35, 'actual_volume_fraction': record['volume_fraction'],
            'angle_degrees': angle, 'load_magnitude_N': 1000., 'mesh_h_mm': record['mesh_h_mm'],
            'status': 'kinematic_diagnostic_only_successful_q4_partial_audit',
            'source_status': record['status'], 'source_q4_solver_passed': True,
            'source_independent_verified': record['independent_verified'],
            'source_independent_failure': record.get('independent_failure'),
            'source_fields_path': relative(fields_path), 'source_fields_sha256': input_hashes[relative(fields_path)],
            'max_displacement_mm': maximum, 'max_pad_displacement_mm': pad_maximum,
            'max_displacement_over_120mm': maximum/120., 'max_pad_displacement_over_120mm': pad_maximum/120.,
            'kinematics': summary, 'factorizations': 0, 'rhs_solves': 0,
            'acceptance_effect': 'None; partial Q4 postprocessing does not confer P2 or reference acceptance'})
        print(json.dumps({'angle': angle, 'maximum_mm': maximum, 'pad_mm': pad_maximum,
            'global_Q_over_e': summary['global']['quadratic_to_linear_area_l2_ratio'],
            'regular_Q_over_e': summary['regular']['quadratic_to_linear_area_l2_ratio'],
            'polar_angle_p95_degrees': np.rad2deg(summary['global']['polar_rotation_absolute_radians']['area_weighted_p95'])}), flush=True)
    for path in input_paths:
        if digest(path) != input_hashes[relative(path)]:
            raise RuntimeError('An input archive changed during read-only postprocessing')
    write_new(output/'manifest.json', {'schema_version': 'prometheus.phase11.repair-kinematics.v1',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'command': sys.argv,
        'catalog_key': entry['catalog_key'], 'catalog_identity_sha256': entry['catalog_identity_sha256'],
        'input_sha256': input_hashes, 'code_protocol_sha256': code_hashes,
        'python': sys.version, 'numpy': np.__version__, 'source_preservation_checked': True,
        'basis': 'RIGHT/DOWN columns already scaled to1000N; angle coefficients sin(theta),cos(theta)',
        'scope': 'Late computed Q4 addendum only; no mechanics solve or acceptance change'})
    write_new(output/'records.json', rows)
    flat = [flatten_summary(row) for row in rows]
    with (output/'summary.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    cost = {'wall_seconds': perf_counter()-started, 'cpu_seconds': process_time()-cpu,
            'factorizations': 0, 'rhs_solves': 0, 'mechanics_solves': 0,
            'records': len(rows), 'source_attempt_wall_seconds': record.get('attempt_wall_seconds'),
            'note': 'Postprocessing cost is separate from the archived Q4/P2 attempt cost; the P2 failure remains recorded'}
    write_new(output/'cost.json', cost)
    print(json.dumps(cost, indent=2), flush=True)
    return rows, cost


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--attempt', default='results/phase11/refinement/attempts/v035_repair/h0.25')
    args = parser.parse_args()
    run(args.output, args.attempt)


if __name__ == '__main__':
    main()
