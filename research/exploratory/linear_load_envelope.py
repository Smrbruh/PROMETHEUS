"""Algebra-only Phase 0 check; not a bracket FE solver or a performance result.

Run: python research/exploratory/linear_load_envelope.py
Writes a deterministic report next to this file. Requires NumPy only.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np


def envelope(q, low, high):
    """Max a(theta)^T Q a(theta), a=[sin(theta), cos(theta)]."""
    # c'=(Q00-Q11) sin(2 theta)+2 Q01 cos(2 theta).
    a, b = q[0, 0] - q[1, 1], 2 * q[0, 1]
    angles = [low, high]
    if abs(a) + abs(b) > 0:
        base = np.arctan2(-b, a)
        kmin = int(np.floor((2 * low - base) / np.pi)) - 1
        kmax = int(np.ceil((2 * high - base) / np.pi)) + 1
        angles += [float((base + k * np.pi) / 2)
                   for k in range(kmin, kmax + 1)
                   if low <= (base + k * np.pi) / 2 <= high]
    values = [np.array([np.sin(t), np.cos(t)]) @ q @
              np.array([np.sin(t), np.cos(t)]) for t in angles]
    winner = int(np.argmax(values))
    return float(values[winner]), angles[winner]


def main():
    rng = np.random.default_rng(1729)
    worst_reconstruction = worst_grid_gap = 0.0
    interior_cases = 0
    for _ in range(40):
        # Deliberately synthetic SPD system: arbitrary units, no geometry.
        a = rng.normal(size=(12, 12))
        k = a.T @ a + 2 * np.eye(12)
        f = rng.normal(size=(12, 2))
        q = f.T @ np.linalg.solve(k, f)
        q = (q + q.T) / 2
        lo, hi = np.deg2rad([-20., 20.])
        exact, argmax = envelope(q, lo, hi)
        interior_cases += int(lo + 1e-12 < argmax < hi - 1e-12)
        theta = np.linspace(lo, hi, 10001)
        coeff = np.array([np.sin(theta), np.cos(theta)])
        reduced = np.einsum('ij,ij->j', coeff, q @ coeff)
        direct_f = f @ coeff[:, ::100]
        direct = np.einsum('ij,ij->j', direct_f,
                           np.linalg.solve(k, direct_f))
        reconstruction = np.max(np.abs(direct - reduced[::100])) / exact
        grid_gap = abs(exact - reduced.max()) / exact
        worst_reconstruction = max(worst_reconstruction, reconstruction)
        worst_grid_gap = max(worst_grid_gap, grid_gap)
        assert reconstruction < 1e-12
        assert grid_gap < 1e-7
        assert exact >= reduced.max() - 1e-12 * exact
        # Uniform material scaling requires all stiffness terms to scale.
        load_scale, modulus_scale = 1.1, .95
        base_force = f[:, 1]
        base_c = base_force @ np.linalg.solve(k, base_force)
        scaled_f = load_scale * base_force
        scaled_c = scaled_f @ np.linalg.solve(modulus_scale * k, scaled_f)
        assert np.isclose(scaled_c / base_c, load_scale**2 / modulus_scale,
                          rtol=1e-12, atol=0)
    area, thickness, density = 120 * 40 + 40 * 40, 5, 2.70e-6
    pad_area = 40 * 10 + 10 * 20
    report = {
        'scope': 'synthetic linear algebra and benchmark bookkeeping only',
        'bracket_fea_or_optimization_run': False,
        'seed': 1729, 'synthetic_systems': 40,
        'dimension_per_system': 12,
        'load_angle_degrees': [-20, 20],
        'max_relative_direct_vs_two_rhs_error': worst_reconstruction,
        'max_relative_exact_vs_dense_grid_gap': worst_grid_gap,
        'interior_angle_maxima_count': interior_cases,
        'uniform_load_1_1_and_modulus_0_95_compliance_ratio': 1.1**2/.95,
        'benchmark_bookkeeping': {
            'area_mm2': area, 'solid_volume_mm3': area*thickness,
            'solid_mass_kg': area*thickness*density,
            'volume_fraction': .4, 'budget_volume_mm3': .4*area*thickness,
            'budget_mass_kg': .4*area*thickness*density,
            'preserved_volume_mm3': pad_area*thickness,
            'variable_budget_volume_mm3': (.4*area-pad_area)*thickness,
            'q4_element_counts_h2_h1_h0_5': [int(area/h**2) for h in [2,1,.5]],
        },
        'all_assertions_passed': True,
        'python': platform.python_version(), 'numpy': np.__version__,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    target = Path(__file__).with_name('linear_load_envelope_results.json')
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
