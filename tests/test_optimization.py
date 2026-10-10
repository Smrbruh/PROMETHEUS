"""Tests of the actual pad-aware filter, derivatives and evaluated-state OC."""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest

from prometheus.geometry import make_l_mesh, bracket_pads, bracket_load_basis, bracket_fixed_dofs, subset_mesh
from prometheus.models import Material, Benchmark, NumericalError
from prometheus import fem
from prometheus.optimization import (DensityFilter, compliance_sensitivity, oc_update,
    solve_simp, simple_l_design, binary_by_rank, geometry_status, sensitivity_audit)


def test_filter_matches_independently_dense_formula_and_pad_override():
    mesh = make_l_mesh(2)
    pads = bracket_pads(mesh)
    filt = DensityFilter(mesh, pads, 6)
    raw = np.linspace(.2, .8, len(mesh.cells))
    raw[pads] = 1
    # Independently calculate selected external-boundary, pad and cutout rows.
    selected = [0, 59, 1200, len(raw)-1]
    centers = mesh.centroids
    for location in ((20, 68), (108, 20), (38, 42), (42, 38)):
        selected.append(int(np.argmin(np.linalg.norm(centers-location, axis=1))))
    physical = filt.physical(raw)
    for i in selected:
        weights = np.maximum(0, 6-np.linalg.norm(centers-centers[i], axis=1))*mesh.areas
        expected = 1 if pads[i] else weights @ raw / weights.sum()
        assert physical[i] == pytest.approx(expected, abs=2e-15)
    np.testing.assert_allclose(np.asarray(filt.H.sum(axis=1)).ravel(), 1, rtol=0, atol=1e-15)
    assert filt.H.shape == (1600, 1600)  # the cutout is absent
    assert np.all(physical[pads] == 1)
    assert not np.allclose(np.asarray(filt.H.sum(axis=0)).ravel(), 1)


def test_initialization_physical_volume_and_exact_chain_rule():
    mesh = make_l_mesh(2)
    pads = bracket_pads(mesh)
    filt = DensityFilter(mesh, pads, 6)
    raw = filt.initialize(.4, .001)
    assert filt.volume(raw, 5) == pytest.approx(12800, abs=1e-8)
    assert np.all(raw[pads] == 1)
    assert np.ptp(raw[~pads]) == 0
    assert raw[~pads][0] != pytest.approx(1960/5800, abs=1e-4)
    gradient = filt.volume_gradient(5)
    assert np.all(gradient[pads] == 0)
    assert np.ptp(gradient[~pads]) > 1  # normalized boundaries cannot use constant h_j
    rng = np.random.default_rng(180)
    direction = rng.normal(size=len(raw))
    direction[pads] = 0
    direction /= np.linalg.norm(direction)
    measured = (filt.volume(raw+.001*direction, 5)-filt.volume(raw-.001*direction, 5))/.002
    assert measured == pytest.approx(gradient @ direction, rel=1e-8, abs=1e-8)
    # Changing only raw pad entries has no physical effect.
    changed = raw.copy()
    changed[pads] = .2
    np.testing.assert_array_equal(filt.physical(changed), filt.physical(raw))
    with pytest.raises(ValueError, match="minimum"):
        filt.initialize(.01)


@pytest.fixture(scope="module")
def derivative_audit():
    return sensitivity_audit(make_l_mesh(2), Material(), Benchmark())


def test_ten_direction_sensitivity_sweeps(derivative_audit):
    audit = derivative_audit
    assert audit["directions"].shape == (10, 1600)
    np.testing.assert_allclose(np.linalg.norm(audit["directions"], axis=1), 1, atol=1e-14)
    assert np.all(audit["directions"][:, bracket_pads(make_l_mesh(2))] == 0)
    assert np.all(audit["compliance_two_adjacent_pass"])
    # Nonlinear central differences must show a decreasing truncation-error
    # region, not merely one fortuitously small step.
    errors = audit["compliance_relative_errors"]
    assert np.all(errors[:, 1] < errors[:, 0])
    # Linear volume has no truncation error; its largest steps minimize roundoff.
    assert np.max(audit["volume_relative_errors"][:, :3]) < 1e-7
    assert audit["factorizations"] == 101


def test_signed_directions_satisfy_adjacent_step_accuracy():
    audit = sensitivity_audit(make_l_mesh(2), Material(), Benchmark(), signed=True)
    assert np.all(audit["compliance_two_adjacent_pass"])
    assert np.all(np.any(audit["directions"][3:] < 0, axis=1))
    assert np.all(np.any(audit["directions"][3:] > 0, axis=1))
    # Sign cancellation makes the derivative small and increases subtraction
    # error; monotonic error decrease is not a property below roundoff. The
    # same 1e-4 accuracy at two adjacent steps is still mandatory here.
    assert np.max(audit["volume_relative_errors"][:, :3]) < 1e-7


def test_oc_feasible_endpoint_move_limit_and_invalid_sensitivities():
    mesh = make_l_mesh(2)
    filt = DensityFilter(mesh, bracket_pads(mesh), 6)
    raw = filt.initialize(.4)
    gradient = -np.linspace(1, 20, len(raw))
    updated = oc_update(raw, gradient, filt, 5)
    assert filt.volume(updated, 5) <= 12800
    assert 12800-filt.volume(updated, 5) < 1e-8*32000
    assert np.max(np.abs(updated-raw)) <= .2+1e-15
    assert np.all(updated[filt.pads] == 1)
    assert np.min(updated) >= .001
    with pytest.raises(ValueError, match="positive sensitivity"):
        oc_update(raw, -gradient, filt, 5)
    with pytest.raises(ValueError, match="finite"):
        oc_update(raw, np.full(len(raw), np.nan), filt, 5)
    with pytest.raises(ValueError, match="move-limit"):
        oc_update(np.ones(len(raw)), gradient, filt, 5)
    zero_update = oc_update(raw, np.zeros(len(raw)), filt, 5)
    assert np.isfinite(zero_update).all()
    assert filt.volume(zero_update, 5) <= 12800


def test_optimizer_budget_returns_last_evaluated_state_and_repeats():
    mesh, material = make_l_mesh(10), Material()
    benchmark = Benchmark(max_evaluations=4)
    a = solve_simp(mesh, material, benchmark)
    b = solve_simp(mesh, material, benchmark)
    assert a.status == b.status == "budget_exhausted"
    assert a.evaluations == 4
    assert len(a.history) == 4
    np.testing.assert_array_equal(a.raw, b.raw)
    np.testing.assert_array_equal(a.physical, b.physical)
    np.testing.assert_array_equal([r["compliance"] for r in a.history],
                                  [r["compliance"] for r in b.history])
    modulus = material.E*(benchmark.emin_ratio+(1-benchmark.emin_ratio)*a.physical**3)
    solved = fem.solve(mesh, modulus, 1000*bracket_load_basis(mesh)[:, 1],
                       bracket_fixed_dofs(mesh), material.nu, material.thickness)
    assert solved.compliance[0] == pytest.approx(a.history[-1]["compliance"], rel=1e-12)
    assert a.history[a.best_evaluation-1]["compliance"] == min(r["compliance"] for r in a.history)
    assert np.all(a.physical[bracket_pads(mesh)] == 1)


def test_optimizer_never_updates_after_single_evaluation(monkeypatch):
    import prometheus.optimization as optimization
    mesh, material = make_l_mesh(10), Material()
    benchmark = Benchmark(max_evaluations=1)
    raw = DensityFilter(mesh, bracket_pads(mesh), 6).initialize(.4)
    def forbidden_update(*args, **kwargs):
        raise AssertionError("An unevaluated final update was attempted")
    monkeypatch.setattr(optimization, "oc_update", forbidden_update)
    result = solve_simp(mesh, material, benchmark)
    np.testing.assert_array_equal(result.raw, raw)
    assert result.evaluations == 1


def test_convergence_requires_minimum_count_and_full_history(monkeypatch):
    import prometheus.optimization as optimization
    mesh, material = make_l_mesh(10), Material()
    monkeypatch.setattr(optimization, "oc_update", lambda raw, *args: raw.copy())
    result = solve_simp(mesh, material, Benchmark(max_evaluations=35))
    assert result.status == "converged"
    assert result.evaluations == 30
    assert len(result.history) == 30
    assert all(row["max_change"] == 0 for row in result.history[-10:])
    assert result.history[-1]["relative_compliance_spread"] == 0


def test_solver_failure_retains_history_and_is_not_a_result(monkeypatch):
    import prometheus.optimization as optimization
    mesh, material = make_l_mesh(10), Material()
    original = fem.solve
    calls = 0
    def failure_after_one(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise NumericalError("injected singular system", code="singular_system")
        return original(*args, **kwargs)
    monkeypatch.setattr(optimization.fem, "solve", failure_after_one)
    with pytest.raises(NumericalError, match="singular system") as failed:
        solve_simp(mesh, material, Benchmark(max_evaluations=4))
    assert len(failed.value.optimization_history) == 1
    assert failed.value.last_evaluated_raw is not None
    assert failed.value.failed_raw is not None


def test_binary_construction_exact_volume_tie_break_and_connectivity():
    mesh = make_l_mesh(2)
    pads = bracket_pads(mesh)
    simple = simple_l_design(mesh)
    assert simple.sum() == 640
    assert simple[pads].sum() == 150
    assert geometry_status(mesh, simple)["valid"]
    equal = binary_by_rank(mesh, np.full(len(mesh.cells), .5))
    free = np.flatnonzero(~pads)
    expected = pads.copy()
    expected[free[:490]] = True
    np.testing.assert_array_equal(equal, expected)
    with pytest.raises(ValueError, match="integer cell"):
        binary_by_rank(mesh, np.ones(len(mesh.cells)), .4001)
    with pytest.raises(ValueError, match="smaller"):
        simple_l_design(mesh, .05)


def test_connectivity_rejects_vertex_join_and_missing_pad():
    mesh = make_l_mesh(2)
    # Two interior cells share only one vertex: they must remain two components.
    binary = np.zeros(len(mesh.cells), dtype=bool)
    binary[0] = True
    binary[61] = True
    status = geometry_status(mesh, binary)
    assert status["component_count"] == 2
    assert not status["connected"]
    assert status["offending_area"] == 8
    simple = simple_l_design(mesh)
    simple[np.flatnonzero(bracket_pads(mesh))[0]] = False
    assert not geometry_status(mesh, simple)["pads_preserved"]
    # Even a mesh with the missing pad removed must not appear complete.
    trimmed = subset_mesh(mesh, simple)
    assert not geometry_status(trimmed, np.ones(len(trimmed.cells)))["valid"]


@pytest.mark.parametrize("radius", [0, -1, np.nan])
def test_filter_invalid_radius(radius):
    mesh = make_l_mesh(10)
    with pytest.raises(ValueError):
        DensityFilter(mesh, bracket_pads(mesh), radius)


def test_invalid_initial_and_filter_input():
    mesh = make_l_mesh(10)
    filt = DensityFilter(mesh, bracket_pads(mesh), 6)
    with pytest.raises(ValueError, match="budget"):
        solve_simp(mesh, Material(), Benchmark(max_evaluations=1), np.ones(len(mesh.cells)))
    with pytest.raises(ValueError, match="bounds"):
        solve_simp(mesh, Material(), Benchmark(max_evaluations=1), np.zeros(len(mesh.cells)))
    with pytest.raises(ValueError, match="finite vector"):
        filt.physical(np.zeros(3))
    with pytest.raises(ValueError, match="Thickness"):
        filt.volume_gradient(0)
