"""Phase 1.1 labels use evidence, not optimistic choice of a finest mesh."""
from copy import deepcopy
import pytest

from prometheus.reference_labels import reference_label


def records(levels=(.5, .25, .125)):
    result = []
    for h, compliance, p2 in zip(levels, (97., 99., 99.5), (98., 99.3, 99.8)):
        result.append({
            "candidate": "frozen_binary", "hash": "same-array-hash", "representation": "binary",
            "problem_signature": "same-physics", "source_geometry_hash": "same-geometry",
            "source_h_mm": 2., "mesh_h_mm": h, "status": "evaluated", "pads_preserved": True,
            "geometry_status": {"valid": True}, "volume_fraction": .4, "volume_budget_fraction": .4,
            "solver_passed": True, "compliance_n_mm": compliance, "max_pad_displacement_mm": compliance/1000,
            "load_gram": [[50., -5.], [-5., compliance]],
            "diagnostics": {"passed": True, "residual_max": 1e-11, "symmetry_error": 1e-16,
                            "energy_error_max": 1e-11, "element_energy_error_max": 1e-11,
                            "force_balance_max": 1e-10, "moment_balance_max": 1e-10,
                            "minimum_relative_lu_pivot": .001},
            "independent": {
                "solver_passed": True, "compliance_n_mm": p2, "max_pad_displacement_mm": p2/1000,
                "load_gram": [[50., -5.], [-5., p2]],
                "q4_relative_compliance_difference": 0., "q4_relative_pad_difference": 0.,
                "diagnostics": {"passed": True, "free_residual_relative": [1e-11, 2e-11],
                                "symmetry_error": 1e-16, "energy_error_relative": [1e-11, 2e-11],
                                "resultant_balance_relative": [1e-10, 2e-10],
                                "moment_balance_relative": [1e-10, 2e-10],
                                "lu_pivot_min_max_ratio": .001}}})
    return result


def geometry():
    return {"valid": True, "edge_connected": True, "pads_preserved": True,
            "local_point_contacts_checked": True, "local_point_contact_count": 0,
            "candidate_hash": "same-array-hash"}


def test_three_new_nested_levels_accept_without_claiming_truth_or_changing_history():
    rows = records()
    old = {"passed": False, "status": "resolution_unverified", "old_levels": [2., 1., .5]}
    original_rows, original_old = deepcopy(rows), deepcopy(old)
    actual = reference_label(rows, historical_label=old, geometry_audit=geometry())
    assert actual["selected_levels_mm"] == [.5, .25, .125]
    assert actual["passed"] and actual["eligible_reference"]
    assert actual["classification"] == "numerically_accepted"
    assert actual["constraint_status"] == "within_budget"
    assert not actual["finest_is_exact_truth"]
    assert actual["physical_validation"] == "not_assessed"
    assert actual["historical_label"] == old == original_old
    assert rows == original_rows
    actual["historical_label"]["passed"] = True
    assert old == original_old


def test_use_final_three_not_coarse_history_to_relax_a_failing_refinement():
    rows = records()
    old = deepcopy(rows[0]); old["mesh_h_mm"] = 1.; old["compliance_n_mm"] = 1.
    rows.insert(0, old)
    rows[1]["compliance_n_mm"] = 98.9
    actual = reference_label(rows, geometry_audit=geometry())
    assert actual["selected_levels_mm"] == [.5, .25, .125]
    assert not actual["passed"]
    assert not actual["gates"]["q4_decreasing_compliance_changes"]


@pytest.mark.parametrize("status", ["failed", "timeout", "resource_cap", "cost_cap", "preflight_skipped", "ill_conditioned"])
def test_latest_failed_or_capped_attempt_cannot_be_replaced_by_older_success(status):
    rows = records((1., .5, .25))
    failure = deepcopy(rows[-1]); failure.update(mesh_h_mm=.125, status=status)
    # Even stale success flags or metrics cannot turn a failed attempt into a solve.
    rows.append(failure)
    label = reference_label(rows, geometry_audit=geometry())
    assert label["selected_levels_mm"] == [.5, .25, .125]
    assert not label["passed"] and not label["eligible_reference"]
    assert "unsuccessful_selected_attempt" in label["reason_codes"]
    older = reference_label(rows, planned_levels=[1., .5, .25], geometry_audit=geometry())
    assert "finer_attempt_omitted" in older["reason_codes"]
    assert not older["passed"]


def test_latest_planned_but_missing_level_is_not_finest_available_truth():
    rows = records((1., .5, .25))
    label = reference_label(rows, planned_levels=[1., .5, .25, .125], historical_label={"passed": True})
    assert label["selected_levels_mm"] == [.5, .25, .125]
    assert "missing_committed_level" in label["reason_codes"]
    assert not label["passed"] and label["reference_estimate"] is None


@pytest.mark.parametrize("levels", [(.5, .125, .0625), (.6, .3, .15), (.5, .25, .25)])
def test_non_factor_two_non_source_nested_or_duplicate_meshes_fail(levels):
    label = reference_label(records(levels), geometry_audit=geometry())
    assert not label["passed"]


@pytest.mark.parametrize("key", ["candidate", "hash", "representation", "problem_signature", "source_geometry_hash", "source_h_mm"])
def test_same_problem_provenance_is_required_across_all_rows(key):
    rows = records()
    rows[0][key] = 1. if key == "source_h_mm" else "different"
    assert not reference_label(rows, geometry_audit=geometry())["passed"]
    rows = records()
    for row in rows:
        del row[key]
    assert not reference_label(rows, geometry_audit=geometry())["passed"]


@pytest.mark.parametrize("value", [float("nan"), float("inf"), 0., -1., None])
@pytest.mark.parametrize("independent", [False, True])
def test_invalid_nominal_metrics_never_accept(value, independent):
    rows = records()
    part = rows[-1]["independent"] if independent else rows[-1]
    part["compliance_n_mm"] = value
    assert not reference_label(rows, geometry_audit=geometry())["passed"]


@pytest.mark.parametrize("field,value", [("residual_max", 1e-7), ("energy_error_max", float("nan")),
                                          ("minimum_relative_lu_pivot", 1e-15)])
def test_q4_stale_pass_flag_cannot_hide_invalid_diagnostics(field, value):
    rows = records(); rows[-1]["diagnostics"][field] = value
    label = reference_label(rows, geometry_audit=geometry())
    assert not label["passed"] and "solver_diagnostics_failed" in label["reason_codes"]


@pytest.mark.parametrize("field,value", [("free_residual_relative", [1e-7]),
                                          ("resultant_balance_relative", []),
                                          ("lu_pivot_min_max_ratio", 1e-15)])
def test_p2_stale_pass_flag_cannot_hide_invalid_diagnostics(field, value):
    rows = records(); rows[-1]["independent"]["diagnostics"][field] = value
    assert not reference_label(rows, geometry_audit=geometry())["passed"]


def test_independent_solution_must_use_same_h_and_exist_at_all_three_levels():
    rows = records(); rows[-1]["independent"]["mesh_h_mm"] = .25
    assert not reference_label(rows, geometry_audit=geometry())["passed"]
    rows = records(); del rows[0]["independent"]
    assert not reference_label(rows, geometry_audit=geometry())["passed"]


def test_finest_agreement_is_recomputed_not_read_from_cached_zero():
    rows = records()
    for row in rows:
        row["independent"]["compliance_n_mm"] *= 2
    label = reference_label(rows, geometry_audit=geometry())
    assert not label["passed"]
    assert not label["gates"]["finest_same_h_compliance_agreement"]


@pytest.mark.parametrize("independent", [False, True])
def test_compliance_must_not_decrease_even_when_absolute_changes_shrink(independent):
    rows = records()
    for row, value in zip(rows, (100., 99.9, 99.85)):
        (row["independent"] if independent else row)["compliance_n_mm"] = value
    label = reference_label(rows, geometry_audit=geometry())
    assert not label["passed"]
    assert not label["gates"]["p2_compliance_monotonicity" if independent else "q4_compliance_monotonicity"]


def test_pad_is_not_required_to_be_monotone_but_final_change_is_gated():
    rows = records()
    for row, value in zip(rows, (.101, .1, .1005)):
        row["max_pad_displacement_mm"] = value
        row["independent"]["max_pad_displacement_mm"] = value
    assert reference_label(rows, geometry_audit=geometry())["passed"]
    rows[-1]["max_pad_displacement_mm"] = .11
    assert not reference_label(rows, geometry_audit=geometry())["passed"]


def test_over_budget_is_separate_from_numerical_acceptance():
    rows = records()
    for row in rows:
        row["volume_fraction"] = .405
        row["volume_feasible"] = True  # obsolete cached claim must not decide eligibility
    label = reference_label(rows, geometry_audit=geometry())
    assert label["passed"] and label["numerical_status"] == "numerically_accepted"
    assert label["constraint_status"] == "over_budget"
    assert label["classification"] == "constraint_over_budget"
    assert not label["eligible_reference"]


def test_disconnected_or_missing_pad_geometry_is_separately_rejected():
    rows = records(); rows[0]["geometry_status"]["valid"] = False
    label = reference_label(rows, geometry_audit=geometry())
    assert label["geometry_status"] == label["classification"] == "geometry_rejected"
    assert not label["eligible_reference"]
    rows = records(); audit = geometry(); audit["pads_preserved"] = False
    assert reference_label(rows, geometry_audit=audit)["classification"] == "geometry_rejected"


def test_globally_connected_local_point_contacts_withhold_eligibility_not_repair_geometry():
    audit = geometry(); audit["local_point_contact_count"] = 3
    label = reference_label(records(), geometry_audit=audit)
    assert label["passed"]  # response of the unchanged tied-node model can converge
    assert label["geometry_status"] == "unresolved_point_contacts"
    assert not label["eligible_reference"]
    assert label["classification"] == "unresolved"
    assert audit["local_point_contact_count"] == 3


def test_missing_or_wrong_identity_contact_audit_never_means_valid_geometry():
    label = reference_label(records())
    assert label["passed"] and not label["eligible_reference"]
    assert label["geometry_status"] == "geometry_unassessed"
    audit = geometry(); audit["candidate_hash"] = "other-candidate"
    assert not reference_label(records(), geometry_audit=audit)["eligible_reference"]


def test_geometry_schema_matches_binary_audit_diagnostics():
    audit = {"geometry_valid": True, "edge_connected": True, "pads_preserved": True,
             "local_point_contact_detection_checked": True, "local_point_contact_count": 0}
    assert reference_label(records(), geometry_audit=audit)["eligible_reference"]
    audit["local_point_contact_count"] = 2
    assert reference_label(records(), geometry_audit=audit)["geometry_status"] == "unresolved_point_contacts"
    audit["geometry_valid"] = False
    assert reference_label(records(), geometry_audit=audit)["geometry_status"] == "geometry_rejected"


def test_load_angle_is_part_of_reference_identity_and_nominal_default_is_explicit():
    rows = records()
    assert reference_label(rows, geometry_audit=geometry())["response_angle_degrees"] == 0.
    for row in rows:
        row["response_angle_degrees"] = -20.
    label = reference_label(rows, geometry_audit=geometry())
    assert label["passed"] and label["response_angle_degrees"] == -20.
    assert label["frozen_identity"]["response_angle_degrees"] == -20.
    rows[-1]["response_angle_degrees"] = 20.
    assert "mixed_response_setting" in reference_label(rows)["reason_codes"]
    rows[-1]["response_angle_degrees"] = float("nan")
    assert "invalid_response_setting" in reference_label(rows)["reason_codes"]


def test_raw_gram_reciprocity_checked_instead_of_inheriting_stale_pass_flag():
    rows = records()
    rows[-1]["independent"]["load_gram"] = [[50., -5.1], [-5., 99.8]]
    label = reference_label(rows, geometry_audit=geometry())
    assert not label["passed"] and "solver_diagnostics_failed" in label["reason_codes"]


@pytest.mark.parametrize("rows", [None, [], ["invalid"], [{"status": "cost_cap", "mesh_h_mm": .125}]])
def test_malformed_or_incomplete_attempts_report_unresolved(rows):
    label = reference_label(rows, historical_label={"passed": True})
    assert not label["passed"]
    assert label["reference_estimate"] is None


def test_original_inequality_is_strict_at_one_percent():
    rows = records()
    for row, value in zip(rows, (96., 99., 100.)):
        row["compliance_n_mm"] = value
        row["independent"]["compliance_n_mm"] = value
    label = reference_label(rows, geometry_audit=geometry())
    assert not label["passed"]
    assert label["q4_refinement"]["last_compliance_change"] == .01


def test_roundoff_sized_compliance_decline_is_allowed():
    rows = records()
    for row in rows:
        row["compliance_n_mm"] = 100.
        row["independent"]["compliance_n_mm"] = 100.
        row["max_pad_displacement_mm"] = .1
        row["independent"]["max_pad_displacement_mm"] = .1
    rows[0]["compliance_n_mm"] += 5e-7
    assert reference_label(rows, geometry_audit=geometry())["passed"]
