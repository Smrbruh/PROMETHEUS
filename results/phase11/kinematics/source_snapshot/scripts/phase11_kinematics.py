"""Postprocess archived Q4 fields without performing a mechanics solve.

Original catalogue candidates are read at h=0.5 mm and angles -20,0,20 degrees.
The two geometry-failed .35 reconstructions retain explicit missing-field
records. Canonical B3 density context is separate and uses only rho>=.9 cells.
Existing evidence is never overwritten; the output directory must not exist.
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
import platform
import sys
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
from prometheus.kinematic_audit import audit_q4_kinematics


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def write_new(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def flatten_summary(record):
    flat = {key: record.get(key) for key in ('catalog_key', 'candidate_name', 'volume_setting',
        'angle_degrees', 'status', 'mesh_h_mm', 'max_displacement_mm', 'max_pad_displacement_mm',
        'max_displacement_over_120mm', 'max_pad_displacement_over_120mm', 'source_fields_path')}
    for region in ('global', 'regular'):
        group = record.get('kinematics', {}).get(region, {})
        flat[f'{region}_area_mm2'] = group.get('area_mm2')
        flat[f'{region}_quadratic_to_linear_l2'] = group.get('quadratic_to_linear_area_l2_ratio')
        for metric in ('linear_strain_frobenius', 'quadratic_correction_frobenius',
                       'displacement_gradient_frobenius', 'skew_gradient_frobenius'):
            for statistic in ('area_l2', 'area_rms', 'area_weighted_p95', 'maximum'):
                flat[f'{region}_{metric}_{statistic}'] = (group.get('metrics') or {}).get(metric, {}).get(statistic)
    return flat


def run(output, catalogue, repair_fields=None):
    output = Path(output).resolve()
    # A preexisting empty directory is also refused, making this append-only.
    output.mkdir(parents=True, exist_ok=False)
    started, cpu = perf_counter(), process_time()
    catalogue = Path(catalogue).resolve()
    inputs = {}

    def track(path):
        path = Path(path).resolve()
        key, value = relative(path), digest(path)
        if key in inputs and inputs[key] != value:
            raise RuntimeError('Source evidence changed during postprocessing: '+key)
        inputs[key] = value
        return value

    track(catalogue)
    catalog = json.loads(catalogue.read_text(encoding='utf-8'))
    protocol = ROOT/'research/REFERENCE_LABEL_PROTOCOL.md'
    sources = [Path(__file__).resolve(), ROOT/'prometheus/kinematic_audit.py',
               ROOT/'prometheus/fem.py', protocol]
    source_hashes = {relative(path): digest(path) for path in sources}
    for path in sources:
        target = output/'source_snapshot'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(path.read_bytes())
    angles = (-20., 0., 20.)
    records, supplemental, omitted = [], [], []

    def postprocess(entry, field_path, history_row, scope='original_binary_catalogue'):
        candidate_started = perf_counter()
        field_hash = track(field_path)
        with np.load(field_path, allow_pickle=False) as archive:
            nodes, cells = archive['nodes'], archive['cells']
            basis, rho = archive['u_basis'], archive['rho']
        if basis.shape != (2*len(nodes), 2) or not np.isfinite(basis).all():
            raise ValueError('Archived displacement basis must have exactly right/down columns')
        if scope == 'original_binary_catalogue' and not np.all(rho == 1):
            raise ValueError('Binary diagnostic requires a material-only mesh')
        if scope == 'density_context':
            cells = cells[rho >= .9]
        if len(cells) == 0:
            raise ValueError('No material cells for the declared kinematic region')
        used = np.unique(cells)
        x, y = nodes.T
        on_pad = (((x >= -1e-10) & (x <= 40+1e-10) & (y >= 70-1e-10) & (y <= 80+1e-10))
                  | ((x >= 110-1e-10) & (x <= 120+1e-10) & (y >= 10-1e-10) & (y <= 30+1e-10)))
        pad_nodes = used[on_pad[used]]
        if not len(pad_nodes):
            raise ValueError('The archived field has no pad nodes')
        result = []
        for angle in angles:
            theta = np.deg2rad(angle)
            u = (np.sin(theta)*basis[:, 0]+np.cos(theta)*basis[:, 1]).reshape(-1, 2)
            kinematics, _ = audit_q4_kinematics(nodes, cells, u)
            magnitudes = np.linalg.norm(u, axis=1)
            maximum, pad_maximum = float(magnitudes[used].max()), float(magnitudes[pad_nodes].max())
            result.append({'catalog_key': entry['catalog_key'], 'candidate_name': entry['candidate_name'],
                'candidate_hash': entry.get('candidate_hash'),
                'catalog_identity_sha256': entry.get('catalog_identity_sha256'),
                'volume_setting': entry['volume_setting'], 'scope': scope,
                'angle_degrees': angle, 'load_magnitude_N': 1000.,
                'basis_convention': 'Archived 1000 N RIGHT and DOWN displacement columns; coefficients sin(theta),cos(theta)',
                'mesh_h_mm': history_row['mesh_h_mm'], 'status': 'kinematic_diagnostic_only',
                'historical_nominal_status': entry.get('historical_label', {}).get('status'),
                'source_fields_path': relative(field_path), 'source_fields_sha256': field_hash,
                'source_solver_passed': history_row.get('solver_passed'),
                'max_displacement_mm': maximum, 'max_pad_displacement_mm': pad_maximum,
                'max_displacement_over_120mm': maximum/120.,
                'max_pad_displacement_over_120mm': pad_maximum/120.,
                'nodal_region': 'retained material cells' if scope != 'density_context' else 'only rho>=0.9 cells; excludes ersatz void nodes',
                'kinematics': kinematics, 'factorizations': 0, 'rhs_solves': 0})
        elapsed = perf_counter()-candidate_started
        for row in result:
            row['three_angle_postprocessing_wall_seconds'] = elapsed
            row['time_accounting'] = 'One shared candidate cost, repeated on its three rows; do not sum these row times'
        return result

    originals = [entry for entry in catalog['entries'] if not entry['is_new_geometry']]
    if len(originals) != 12:
        raise ValueError('Expected exactly twelve original binary catalogue entries')
    for entry in originals:
        for key, hash_key in (('array_path', 'source_input_sha256'),
                              ('history_path', 'historical_history_sha256'),
                              ('source_mesh_path', 'source_mesh_file_sha256')):
            observed = track(ROOT/entry[key])
            if entry.get(hash_key) != observed:
                raise ValueError('Catalogue provenance mismatch: '+entry['catalog_key']+' '+key)
        manifest_path = ROOT/entry['original_run']/'manifest.json'
        track(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest['benchmark']['load'] != 1000:
            raise ValueError('Unexpected archived load magnitude')
        history_path = ROOT/entry['history_path']
        history = json.loads(history_path.read_text(encoding='utf-8'))
        candidates = [row for row in history if row.get('mesh_h_mm') == .5 and row.get('solver_passed')]
        if not entry['geometry']['geometry_valid'] or not candidates:
            status = ('not_evaluated_geometry_failed_no_field' if not entry['geometry']['geometry_valid']
                      else 'not_evaluated_no_valid_h0_5_field')
            for angle in angles:
                records.append({'catalog_key': entry['catalog_key'], 'candidate_name': entry['candidate_name'],
                    'candidate_hash': entry['candidate_hash'], 'volume_setting': entry['volume_setting'],
                    'angle_degrees': angle, 'scope': 'original_binary_catalogue', 'status': status,
                    'historical_nominal_status': entry['historical_label']['status'],
                    'kinematics': None, 'max_displacement_mm': None, 'max_pad_displacement_mm': None,
                    'factorizations': 0, 'rhs_solves': 0,
                    'reason': 'No displacement field is invented for rejected or unevaluated geometry'})
            continue
        row = candidates[0]
        if row.get('hash') != entry['candidate_hash']:
            raise ValueError('Candidate displacement history has a different frozen identity')
        field_path = history_path.parent/'h0.5_fields.npz'
        records.extend(postprocess(entry, field_path, row))
        print('Postprocessed '+entry['catalog_key']+' at three angles', flush=True)

    # Supplemental canonical relaxed model is deliberately outside the binary catalogue.
    density_dir = ROOT/'results/phase1_v040_verified/candidates/B3'
    density_history_path = density_dir/'evaluations.json'
    track(density_history_path)
    density_history = json.loads(density_history_path.read_text(encoding='utf-8'))
    density_row = next(row for row in density_history if row['mesh_h_mm'] == .5)
    density_entry = {'catalog_key': 'v040_B3_density_context', 'candidate_name': 'B3',
                     'candidate_hash': density_row['hash'], 'volume_setting': .4}
    supplemental.extend(postprocess(density_entry, density_dir/'h0.5_fields.npz', density_row, 'density_context'))

    repairs = [entry for entry in catalog['entries'] if entry['is_new_geometry']]
    for entry in repairs:
        # Explicit path selection avoids accidentally adopting a different repair
        # or one of the old parent fields as a new geometry's solution.
        if repair_fields is None:
            omitted.append({'catalog_key': entry['catalog_key'], 'status': 'not_evaluated_no_repair_fields_supplied',
                            'reason': 'Repair audit was unavailable for this independent postprocessing batch; no waiting or parent-field substitution'})
            continue
        field_path = Path(repair_fields).resolve()
        history_path = field_path.with_name(field_path.name.replace('_fields.npz', '.json'))
        track(history_path)
        row = json.loads(history_path.read_text(encoding='utf-8'))
        if row.get('hash') != entry['candidate_hash'] or not row.get('solver_passed'):
            raise ValueError('Supplied repair fields must have matching hash and successful Q4 solve metadata')
        supplemental.extend(postprocess(entry, field_path, row, 'repaired_geometry'))

    for key, value in inputs.items():
        if digest(ROOT/key) != value:
            raise RuntimeError('Read-only evidence changed during postprocessing: '+key)
    manifest = {'schema_version': 'prometheus.phase11.kinematics.v1',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'command': sys.argv,
        'catalogue': relative(catalogue), 'input_sha256': inputs, 'code_protocol_sha256': source_hashes,
        'python': sys.version, 'numpy': np.__version__, 'platform': platform.platform(),
        'angles_degrees': list(angles), 'original_candidate_count': len(originals),
        'source_preservation_checked': True, 'no_mechanics_solves': True,
        'load_normalization': 'No rescaling: supplied archived bases already include 1000 N',
        'regular_region': '6 mm from the five original singular points; new pixel corners are not excluded',
        'acceptance': 'No kinematic threshold, nonlinear equilibrium, safety, or response-error certification'}
    write_new(output/'manifest.json', manifest)
    write_new(output/'records.json', records)
    write_new(output/'supplemental_context.json', supplemental)
    write_new(output/'omitted.json', omitted)
    flat = [flatten_summary(row) for row in records]
    with (output/'summary.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    cost = {'wall_seconds': perf_counter()-started, 'cpu_seconds': process_time()-cpu,
            'factorizations': 0, 'rhs_solves': 0, 'mechanics_solves': 0,
            'original_records': len(records),
            'original_evaluated_records': sum(row['status'] == 'kinematic_diagnostic_only' for row in records),
            'original_unavailable_records': sum(row['status'] != 'kinematic_diagnostic_only' for row in records),
            'supplemental_records': len(supplemental),
            'note': 'Includes file reading, hashing, postprocessing and preceding JSON/CSV writes; excludes earlier archived solve costs'}
    write_new(output/'cost.json', cost)
    print(json.dumps(cost, indent=2), flush=True)
    return records, supplemental, cost


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--catalogue', default='results/phase11/geometry/catalog.json')
    parser.add_argument('--repair-fields', help='Optional already completed repair *_fields.npz with matching adjacent JSON metadata')
    args = parser.parse_args()
    run(args.output, args.catalogue, args.repair_fields)


if __name__ == '__main__':
    main()
