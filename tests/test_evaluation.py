"""Frozen-design semantics and guards against unsupported acceptance claims."""
from copy import deepcopy
import numpy as np
import pytest

import prometheus.evaluation as evaluation
from prometheus.evaluation import (Candidate, transfer_frozen, evaluate_candidate,
                                    verification_summary, solved_load_envelope)
from prometheus.envelope import compliance_envelope
from prometheus.geometry import make_l_mesh, make_rect_mesh, bracket_pads
from prometheus.models import Material, Benchmark, NumericalError
from prometheus.optimization import simple_l_design


def test_frozen_transfer_preserves_each_parent_and_material_volume():
    source, target = make_l_mesh(2), make_l_mesh(1)
    values = np.arange(len(source.cells), dtype=float)/len(source.cells)
    refined = transfer_frozen(source, values, target)
    lookup = {tuple(index): values[i] for i, index in enumerate(source.cell_indices)}
    expected = np.array([lookup[tuple(index//2)] for index in target.cell_indices])
    np.testing.assert_array_equal(refined, expected)
    assert source.areas @ values == pytest.approx(target.areas @ refined, rel=1e-14)
    assert not np.shares_memory(refined, values)
    binary = simple_l_design(source)
    binary_refined = transfer_frozen(source, binary, target)
    assert np.isin(binary_refined, [0, 1]).all()
    assert source.areas @ binary == target.areas @ binary_refined
    np.testing.assert_array_equal(transfer_frozen(source, values, source), values)


def test_frozen_transfer_rejects_coarsening_shape_mismatch_and_outside_domain():
    source = make_l_mesh(2)
    values = np.ones(len(source.cells))
    with pytest.raises(ValueError, match="subdivision"):
        transfer_frozen(source, values, make_l_mesh(4))
    with pytest.raises(ValueError, match="shape"):
        transfer_frozen(source, values[:-1], make_l_mesh(1))
    with pytest.raises(ValueError, match="outside"):
        transfer_frozen(source, values, make_rect_mesh(120, 80, 1))


@pytest.mark.parametrize("values,representation", [
    ([np.nan], "density"), ([np.inf], "density"), ([-.1], "density"),
    ([1.1], "density"), ([[.5]], "density"), ([.5], "binary"), ([1.], "unknown"),
])
def test_candidate_rejects_invalid_fields(values, representation):
    with pytest.raises(ValueError):
        Candidate("invalid", values, representation)


def test_candidate_copies_input_and_hash_excludes_display_name():
    original = np.array([1., .4, .2])
    candidate = Candidate("display name", original, "density")
    original[0] = 0.
    assert candidate.physical[0] == 1.
    with pytest.raises(ValueError):
        candidate.physical[0] = 0.
    assert candidate.content_hash == Candidate("renamed", [1., .4, .2], "density").content_hash
    assert Candidate("same values", [0., 1.], "density").content_hash != Candidate("same values", [0., 1.], "binary").content_hash


@pytest.mark.parametrize("omit_load_pad", [False, True])
def test_failed_binary_geometry_never_calls_a_solver(monkeypatch, omit_load_pad):
    mesh = make_l_mesh(2)
    pixels = bracket_pads(mesh).astype(float)
    if omit_load_pad:
        pixels[mesh.centroids[:, 0] >= 110] = 0.
    def forbidden_solver(*args, **kwargs):
        pytest.fail("Geometry rejected by edge connectivity must not reach a solver")
    monkeypatch.setattr(evaluation, "solve", forbidden_solver)
    monkeypatch.setattr(evaluation, "bracket_verify", forbidden_solver)
    record, arrays = evaluate_candidate(Candidate("disconnected", pixels, "binary"),
                                        mesh, Material(), Benchmark(), independent=True)
    assert record["status"] == "geometry_failed"
    assert record["geometry_status"]["valid"] is False
    assert record["solver_passed"] is False
    assert record["optimization_feasible"] is False
    assert record["independent_verified"] is False
    assert record["factorizations"] == 0 and record["rhs_solves"] == 0
    assert "compliance_n_mm" not in record
    assert arrays == {}


def test_binary_material_mesh_and_volume_constraint_are_distinct_from_verification():
    mesh = make_l_mesh(2)
    candidate = Candidate("B1", simple_l_design(mesh), "binary")
    record, arrays = evaluate_candidate(candidate, mesh, Material(), Benchmark(), include_stress=False)
    assert record["volume_fraction"] == pytest.approx(.4)
    assert record["q4_elements"] == 640
    assert np.all(arrays["rho"] == 1.)
    assert len(np.unique(arrays["cells"])) == len(arrays["nodes"])
    assert record["solver_passed"] and record["optimization_feasible"]
    assert record["max_material_displacement_mm"] == record["max_displacement_all_nodes_mm"]
    assert record["independent_verified"] is False
    assert record["physical_validation"] == "not_assessed"
    assert record["response_feasibility"] == "not_defined_no_response_limit"


def test_full_solid_is_stiffness_reference_and_not_equal_mass_feasible():
    mesh = make_l_mesh(10)
    candidate = Candidate("B0", np.ones(len(mesh.cells)), "density")
    record, _ = evaluate_candidate(candidate, mesh, Material(), Benchmark(), h=10, include_stress=False)
    assert record["solver_passed"]
    assert record["volume_fraction"] == 1.
    assert not record["volume_feasible"] and not record["optimization_feasible"]
    assert not record["independent_verified"]
    assert record["factorizations"] == 1 and record["rhs_solves"] == 2
    assert record["worst_compliance_n_mm"] >= record["compliance_n_mm"]


def test_single_mesh_independent_solver_agreement_is_not_resolution_verification():
    mesh = make_l_mesh(10)
    candidate = Candidate("B0", np.ones(len(mesh.cells)), "density")
    record, arrays = evaluate_candidate(candidate, mesh, Material(), Benchmark(), h=10,
                                        independent=True, include_stress=False)
    assert record["solver_passed"] and record["independent"]["solver_passed"]
    assert record["independent"]["element"] == "P2_triangle"
    assert record["factorizations"] == 2 and record["rhs_solves"] == 4
    assert "p2_nodes" in arrays and "p2_u_basis" in arrays
    for checked in (record,record["independent"]):
        assert checked["load_gram_diagnostics"]["reciprocity_passed"]
        raw=np.asarray(checked["load_gram"])
        measured=np.linalg.norm(raw-raw.T)/np.linalg.norm(raw)
        assert checked["load_gram_diagnostics"]["reciprocity_relative"] == pytest.approx(measured)
    assert not record["independent_verified"]
    assert not verification_summary([record])["passed"]


def test_failed_solve_propagates_instead_of_returning_a_success(monkeypatch):
    mesh = make_l_mesh(10)
    def failed(*args, **kwargs):
        raise NumericalError("Injected singular solve", "singular_system")
    monkeypatch.setattr(evaluation, "solve", failed)
    with pytest.raises(NumericalError, match="Injected"):
        evaluate_candidate(Candidate("failure", np.ones(len(mesh.cells)), "density"),
                            mesh, Material(), Benchmark(), h=10)


def _converged_records():
    rows = []
    for h, compliance, p2 in ((2., 97., 98.), (1., 99., 99.3), (.5, 99.5, 99.8)):
        rows.append({"candidate": "fixture", "hash": "same-frozen-field",
                     "problem_signature": "same-physics", "source_geometry_hash": "same-geometry",
                     "representation": "density", "source_h_mm": 2., "mesh_h_mm": h,
                     "solver_passed": True, "pads_preserved": True,
                     "volume_feasible": True, "optimization_feasible": True,
                     "status": "evaluated", "compliance_n_mm": compliance,
                     "max_pad_displacement_mm": compliance/1000,
                     "independent": {"solver_passed": True, "compliance_n_mm": p2,
                                     "max_pad_displacement_mm": p2/1000,
                                     "q4_relative_compliance_difference": abs(compliance-p2)/p2,
                                     "q4_relative_pad_difference": abs(compliance-p2)/p2}})
    return rows


def test_verification_requires_complete_mesh_sequence_and_both_solvers():
    records = _converged_records()
    assert verification_summary(records)["passed"]
    assert not verification_summary(records[:-1])["passed"]
    records[1]["independent"]["solver_passed"] = False
    assert not verification_summary(records)["passed"]
    records = _converged_records()
    records[1]["solver_passed"] = False
    assert not verification_summary(records)["passed"]


@pytest.mark.parametrize("mutation", ["hash", "representation", "source_h", "duplicate", "pads",
                                      "problem_signature", "source_geometry_hash"])
def test_verification_cannot_combine_different_or_invalid_frozen_designs(mutation):
    records = _converged_records()
    if mutation == "hash":
        records[1]["hash"] = "different-field"
    elif mutation == "representation":
        records[1]["representation"] = "binary"
    elif mutation == "source_h":
        records[1]["source_h_mm"] = 1.
    elif mutation == "duplicate":
        records.append(deepcopy(records[-1]))
    elif mutation == "pads":
        records[1]["pads_preserved"] = False
    else:
        records[1][mutation] = "different-model-or-geometry"
    assert not verification_summary(records)["passed"]


@pytest.mark.parametrize("bad_metric", ["coarse_fine_difference", "solver_disagreement", "worsening_trend"])
def test_verification_rejects_resolution_or_solver_disagreement(bad_metric):
    records = _converged_records()
    if bad_metric == "coarse_fine_difference":
        records[-1]["compliance_n_mm"] = 90.
    elif bad_metric == "solver_disagreement":
        # Both mesh sequences still converge, but to disagreeing values.
        for record in records:
            record["independent"]["compliance_n_mm"] *= 2
    else:
        records[0]["compliance_n_mm"] = 98.9
    assert not verification_summary(records)["passed"]


def test_numerical_verification_does_not_assert_volume_or_physical_feasibility():
    records = _converged_records()
    for record in records:
        record["volume_feasible"] = False
        record["optimization_feasible"] = False
    checked = verification_summary(records)
    assert checked["passed"]  # B0 is allowed as a numerical reference only.
    assert checked["physical_validation"] == "not_assessed"


@pytest.mark.parametrize("bad_value", [-1., 0., np.inf, np.nan])
def test_verification_rejects_nonphysical_or_nonfinite_metrics(bad_value):
    records = _converged_records()
    records[1]["compliance_n_mm"] = bad_value
    assert not verification_summary(records)["passed"]


@pytest.mark.parametrize("metric,cached_key,output_key", [
    ("compliance_n_mm", "q4_relative_compliance_difference", "finest_compliance_agreement"),
    ("max_pad_displacement_mm", "q4_relative_pad_difference", "finest_pad_agreement"),
])
@pytest.mark.parametrize("cached_gap", [0., -.5, 999., np.nan])
def test_solver_agreement_is_recomputed_from_metrics(metric, cached_key, output_key, cached_gap):
    records = _converged_records()
    records[-1]["independent"][cached_key] = cached_gap
    good = verification_summary(records)
    actual = abs(records[-1][metric]-records[-1]["independent"][metric])/records[-1]["independent"][metric]
    assert good["passed"]
    assert good[output_key] == pytest.approx(actual)
    for record in records:
        record["independent"][metric] *= 2.
    bad = verification_summary(records)
    actual = abs(records[-1][metric]-records[-1]["independent"][metric])/records[-1]["independent"][metric]
    assert not bad["passed"]
    assert bad[output_key] == pytest.approx(actual)


@pytest.mark.parametrize("missing_key", ["hash", "representation", "source_h_mm",
                                        "problem_signature", "source_geometry_hash"])
def test_matching_missing_provenance_cannot_establish_a_frozen_design(missing_key):
    records = _converged_records()
    for record in records:
        del record[missing_key]
    assert not verification_summary(records)["passed"]


def test_solved_gram_adapter_preserves_quadratic_response_and_raw_evidence():
    symmetric=np.array([[2., .3], [.3, 3.]])
    gram=symmetric+np.array([[0., 1e-10], [-1e-10, 0.]])
    original=gram.copy()
    interval=np.deg2rad([-20.,20.])
    with pytest.raises(ValueError,match="symmetric"):
        compliance_envelope(gram,*interval)
    result,diagnostics=solved_load_envelope(gram,*interval)
    expected=compliance_envelope(symmetric,*interval)
    assert result.maximum == pytest.approx(expected.maximum,rel=1e-14)
    np.testing.assert_array_equal(gram,original)
    coefficients=np.array([np.sin(np.linspace(*interval,1001)),np.cos(np.linspace(*interval,1001))])
    raw_values=np.einsum("ij,ij->j",coefficients,gram@coefficients)
    symmetric_values=np.einsum("ij,ij->j",coefficients,symmetric@coefficients)
    np.testing.assert_allclose(raw_values,symmetric_values,rtol=1e-14)
    assert diagnostics["reciprocity_relative"] > 1e-12
    assert diagnostics["reciprocity_relative"] < diagnostics["reciprocity_tolerance"]
    assert diagnostics["uses_symmetric_quadratic_part"]


def test_solved_gram_adapter_rejects_nonacceptable_reciprocity():
    gram=np.array([[2.,.3001],[.2999,3.]])
    with pytest.raises(NumericalError,match="reciprocity") as failure:
        solved_load_envelope(gram,-.2,.2)
    assert failure.value.code == "reciprocity_failed"


def test_solved_gram_adapter_does_not_hide_nonpositive_quadratic_energy():
    with pytest.raises(ValueError,match="positive semidefinite"):
        solved_load_envelope(-np.eye(2),-.2,.2)


@pytest.mark.parametrize("solver", ["q4", "p2"])
def test_nested_compliance_decline_is_rejected_even_with_small_refinement_changes(solver):
    records=_converged_records()
    selected=records if solver=="q4" else [r["independent"] for r in records]
    for record,value in zip(selected,[100.,99.9,99.85]):
        record["compliance_n_mm"]=value
    checked=verification_summary(records)
    assert checked["compliance_refinement_trend_improves"]
    assert not checked[f"{solver}_compliance_nondecreasing"]
    assert not checked["passed"]


def test_monotonicity_allows_only_solve_roundoff_and_does_not_require_monotone_pad():
    records=_converged_records()
    for record in records:
        record["compliance_n_mm"]=100.
        record["independent"]["compliance_n_mm"]=100.
    records[0]["compliance_n_mm"]+=5e-7
    for record,displacement in zip(records,[.101,.1005,.100]):
        record["max_pad_displacement_mm"]=displacement
        record["independent"]["max_pad_displacement_mm"]=displacement
    checked=verification_summary(records)
    assert checked["q4_compliance_nondecreasing"] and checked["p2_compliance_nondecreasing"]
    assert checked["passed"]


def test_independent_failure_preserves_completed_primary_cost_and_fields(monkeypatch):
    mesh=make_l_mesh(10)
    def fail_independent(*args,**kwargs):
        raise NumericalError("Injected P2 failure","singular_system")
    monkeypatch.setattr(evaluation,"bracket_verify",fail_independent)
    with pytest.raises(NumericalError,match="Injected P2") as failure:
        evaluate_candidate(Candidate("partial",np.ones(len(mesh.cells)),"density"),
                           mesh,Material(),Benchmark(),h=10,independent=True,include_stress=False)
    record=failure.value.partial_evaluation
    arrays=failure.value.partial_arrays
    assert record["status"]=="independent_verification_failed"
    assert record["solver_passed"] and not record["independent_verified"]
    assert record["factorizations"]==1 and record["rhs_solves"]==2
    assert record["wall_seconds"]>=record["primary_analysis_seconds"]
    # The Windows process CPU clock can quantize a short solve to zero.
    assert np.isfinite(record["cpu_seconds"]) and record["cpu_seconds"]>=0
    assert record["compliance_n_mm"]>0
    assert arrays["u_basis"].shape==(mesh.ndof,2)


def test_independent_postprocessing_failure_charges_both_completed_solves(monkeypatch):
    mesh=make_l_mesh(10)
    real_adapter=evaluation.solved_load_envelope
    calls=0
    def fail_second(*args,**kwargs):
        nonlocal calls
        calls+=1
        if calls==2:
            raise NumericalError("Injected reciprocity failure","reciprocity_failed")
        return real_adapter(*args,**kwargs)
    monkeypatch.setattr(evaluation,"solved_load_envelope",fail_second)
    with pytest.raises(NumericalError,match="reciprocity") as failure:
        evaluate_candidate(Candidate("partial",np.ones(len(mesh.cells)),"density"),
                           mesh,Material(),Benchmark(),h=10,independent=True,include_stress=False)
    record=failure.value.partial_evaluation
    assert record["factorizations"]==2 and record["rhs_solves"]==4
    assert "load_gram" in record["independent"]
    assert not record["independent"]["solver_passed"]
    assert "p2_u_basis" in failure.value.partial_arrays


@pytest.mark.parametrize("failed_function", ["solved_load_envelope", "stress_metrics"])
def test_primary_postprocessing_failure_preserves_completed_solve(monkeypatch,failed_function):
    mesh=make_l_mesh(10)
    def fail_processing(*args,**kwargs):
        raise NumericalError("Injected primary postprocessing failure","postprocessing_failure")
    def forbidden_independent(*args,**kwargs):
        pytest.fail("An earlier postprocessing failure must not reach the independent solve")
    monkeypatch.setattr(evaluation,failed_function,fail_processing)
    monkeypatch.setattr(evaluation,"bracket_verify",forbidden_independent)
    with pytest.raises(NumericalError,match="primary postprocessing") as failure:
        evaluate_candidate(Candidate("partial",np.ones(len(mesh.cells)),"density"),
                           mesh,Material(),Benchmark(),h=10,independent=True,include_stress=True)
    record=failure.value.partial_evaluation
    arrays=failure.value.partial_arrays
    assert record["status"]=="primary_postprocessing_failed"
    assert record["solver_passed"] and not record["independent_verified"]
    assert record["factorizations"]==1 and record["rhs_solves"]==2
    assert record["compliance_n_mm"]>0 and record["diagnostics"]["passed"]
    assert record["wall_seconds"]>=record["primary_analysis_seconds"]
    # Zero is valid below the process-clock resolution; elapsed time is finite.
    assert np.isfinite(record["cpu_seconds"]) and record["cpu_seconds"]>=0
    np.testing.assert_allclose(record["load_gram"],arrays["forces_basis"].T@arrays["u_basis"])
    assert arrays["u_basis"].shape==(mesh.ndof,2)
    assert np.isfinite(arrays["u_basis"]).all()
    assert "p2_u_basis" not in arrays
