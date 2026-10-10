"""Projection algebra and allocation-free topology predictions."""
from copy import deepcopy

import numpy as np
import pytest

from prometheus.audit_responses import project_response, resource_estimate
from prometheus.geometry import make_rect_mesh, subset_mesh
from prometheus.verification import p2_mesh_from_quads


def fixture():
    nodes = np.array([[0., 80.], [120., 20.], [60., 20.], [115., 10.]])
    basis = np.arange(16, dtype=float).reshape(8, 2)/10
    gram = [[3., .25], [.25, 5.]]
    diagnostic = {"passed": True, "retained_value": 12.}
    record = {"status": "evaluated", "problem_signature": "original-model", "solver_passed": True,
              "independent_verified": True, "load_gram": gram, "diagnostics": diagnostic,
              "compliance_n_mm": 5., "max_pad_displacement_mm": 3.,
              "independent": {"solver_passed": True, "load_gram": [[4., .3], [.3, 6.]],
                              "diagnostics": diagnostic, "compliance_n_mm": 6.,
                              "max_pad_displacement_mm": 4.}}
    arrays = {"nodes": nodes, "u_basis": basis,
              "p2_nodes": nodes.copy(), "p2_u_basis": basis*1.1}
    return record, arrays


@pytest.mark.parametrize("angle", [-20., 0., 20., 13.])
def test_response_projection_matches_direct_basis_displacement_and_quadratic(angle):
    record, arrays = fixture()
    original, original_arrays = deepcopy(record), {k: v.copy() for k, v in arrays.items()}
    projected = project_response(record, arrays, angle)
    coefficients = [np.sin(np.deg2rad(angle)), np.cos(np.deg2rad(angle))]
    for prefix, part in (("", projected), ("p2_", projected["independent"])):
        direct = (arrays[prefix+"u_basis"]@coefficients).reshape(-1, 2)
        expected = np.linalg.norm(direct[[0, 1, 3]], axis=1).max()
        assert part["max_pad_displacement_mm"] == pytest.approx(expected, rel=1e-14)
        assert part["compliance_n_mm"] == pytest.approx(np.array(coefficients)@np.asarray(part["load_gram"])@coefficients, rel=1e-14)
    assert projected["response_projection_status"] == "projected_both_solvers"
    assert projected["response_projection_factorizations"] == projected["response_projection_rhs_solves"] == 0
    assert not projected["independent_verified"]  # must audit the actual new response
    assert record == original
    for key in arrays:
        np.testing.assert_array_equal(arrays[key], original_arrays[key])
    assert projected["load_gram"] == record["load_gram"]
    assert projected["diagnostics"] == record["diagnostics"]
    projected["diagnostics"]["retained_value"] = 99
    assert record["diagnostics"]["retained_value"] == 12


def test_response_identity_preserves_nominal_and_deterministically_distinguishes_angles():
    record, arrays = fixture()
    nominal = project_response(record, arrays, 0.)
    negative = project_response(record, arrays, -20.)
    positive = project_response(record, arrays, 20.)
    assert nominal["problem_signature"] == record["problem_signature"]
    assert negative["base_problem_signature"] == record["problem_signature"]
    assert positive["problem_signature"] != negative["problem_signature"] != record["problem_signature"]
    assert project_response(negative, arrays, 20.)["problem_signature"] == positive["problem_signature"]
    assert project_response(positive, arrays, 0.)["problem_signature"] == record["problem_signature"]


@pytest.mark.parametrize("status", ["failed", "timeout", "preflight_skipped", "geometry_failed"])
def test_failure_or_resource_stop_has_no_invented_projected_metrics(status):
    record, arrays = fixture(); record["status"] = status
    projected = project_response(record, None, 20.)
    assert projected["status"] == status
    assert projected["response_projection_status"] == "not_projected_unsuccessful_attempt"
    for part in (projected, projected["independent"]):
        assert "compliance_n_mm" not in part
        assert "max_pad_displacement_mm" not in part
        assert "compliance_n_mm" in part["base_response_metrics"]


def test_missing_fields_fail_closed_and_preserve_raw_metrics_as_history():
    record, arrays = fixture()
    projected = project_response(record, {}, 20.)
    assert projected["status"] == "response_projection_failed"
    assert "compliance_n_mm" not in projected
    assert projected["base_response_metrics"]["compliance_n_mm"] == 5.
    del arrays["p2_u_basis"]
    projected = project_response(record, arrays, 20.)
    assert projected["status"] == "response_projection_failed"
    assert projected["compliance_n_mm"] > 0  # completed primary projection remains evidence
    assert "compliance_n_mm" not in projected["independent"]
    assert not projected["independent"]["solver_passed"]


@pytest.mark.parametrize("angle", [float("nan"), float("inf"), True, None])
def test_invalid_projection_angle_rejected(angle):
    record, arrays = fixture()
    with pytest.raises(ValueError):
        project_response(record, arrays, angle)


def test_malformed_displacements_are_not_reshaped_into_an_unrelated_field():
    record, arrays = fixture(); arrays["u_basis"] = arrays["u_basis"].T
    with pytest.raises(ValueError, match="basis"):
        project_response(record, arrays, 20.)


@pytest.mark.parametrize("factor", [1, 2, 4])
@pytest.mark.parametrize("selection", ["full", "l", "diagonal", "hole"])
def test_topology_prediction_equals_explicit_small_p2_mesh(factor, selection):
    coarse = make_rect_mesh(3., 3., 1.)
    ix, iy = coarse.cell_indices.T
    if selection == "full":
        mask = np.ones(len(ix), dtype=bool)
    elif selection == "l":
        mask = (ix == 0) | (iy == 0)
    elif selection == "diagonal":
        mask = ix == iy
    else:
        mask = (ix != 1) | (iy != 1)
    estimate = resource_estimate(coarse, mask, 1/factor)
    fine = make_rect_mesh(3., 3., 1/factor)
    parent = fine.cell_indices//factor
    lookup = {tuple(ij): active for ij, active in zip(coarse.cell_indices, mask)}
    active = np.array([lookup[tuple(ij)] for ij in parent])
    material = subset_mesh(fine, active)
    p2 = p2_mesh_from_quads(material)
    assert estimate["q4_cells"] == len(material.cells)
    assert estimate["q4_vertices"] == len(material.nodes)
    assert estimate["p2_dofs"] == p2.ndof
    assert estimate["p2_triangles"] == len(p2.cells)
    expected_bytes = len(p2.cells)*144*40+p2.ndof*64
    assert estimate["estimated_p2_coo_bytes"] == expected_bytes
    assert estimate["allowed"]
    assert not estimate["matrix_assembled"]


def test_resource_caps_and_empty_geometry_do_not_claim_a_failed_mechanical_solution():
    mesh = make_rect_mesh(1., 1., 1.)
    huge = resource_estimate(mesh, [1], 1/1024)
    assert not huge["allowed"] and len(huge["reasons"]) == 2
    assert huge["p2_dofs"] > 1_500_000
    assert huge["estimated_p2_coo_bytes"] > 4*1024**3
    assert huge["factorization_memory_estimate_bytes"] is None
    empty = resource_estimate(mesh, [0], .5)
    assert not empty["allowed"] and empty["q4_cells"] == empty["p2_dofs"] == 0


@pytest.mark.parametrize("spacing", [2., .3, 0., -1., float("nan"), True])
def test_non_nested_or_invalid_resource_spacing_rejected(spacing):
    with pytest.raises(ValueError):
        resource_estimate(make_rect_mesh(1., 1., 1.), [1], spacing)


def test_nonbinary_resource_mask_rejected():
    with pytest.raises(ValueError, match="binary"):
        resource_estimate(make_rect_mesh(1., 1., 1.), [.4], .5)
