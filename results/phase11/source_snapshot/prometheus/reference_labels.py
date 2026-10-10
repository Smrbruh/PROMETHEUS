"""Additive Phase 1.1 labels from a rolling three-level frozen-design audit.

No mechanics is executed here, no archived label is modified, and the finest
available value is never presumed to be truth. Numerical acceptance concerns
only declared-angle compliance and nodal maximum pad-displacement magnitude. Geometry,
budget eligibility and physical validation remain separate statements.

``planned_levels`` means the levels committed to the current audit, including
an attempted/skipped terminal level. A future conditional refinement need not
be committed before deciding whether the current three-level audit suffices.
Once a latest level is committed, a failure, missing record or resource stop at
that level cannot be hidden by falling back to three older successful solves.
"""
from copy import deepcopy
import math
from numbers import Real


_IDENTITY = ("candidate", "hash", "representation", "source_h_mm",
             "problem_signature", "source_geometry_hash")
_STRING_IDENTITY = tuple(k for k in _IDENTITY if k != "source_h_mm")
_METRICS = ("compliance_n_mm", "max_pad_displacement_mm")
_MONOTONICITY_TOLERANCE = 1e-8
_TREND_TOLERANCE = 1e-10


def _finite(value, positive=False):
    return (isinstance(value, Real) and not isinstance(value, bool)
            and math.isfinite(value) and (not positive or value > 0))


def _maximum(values):
    """Accept a scalar or nonempty one-dimensional diagnostic sequence."""
    if _finite(values):
        return float(values)
    if isinstance(values, (str, bytes, dict)):
        return None
    try:
        items = list(values)
    except TypeError:
        return None
    if not items or not all(_finite(x) and x >= 0 for x in items):
        return None
    return max(float(x) for x in items)


def _diagnostic_issues(row, independent=False):
    label = "P2" if independent else "Q4"
    issues = []
    if row.get("solver_passed") is not True:
        issues.append(f"{label} solver_passed is not explicitly true")
    diagnostic = row.get("diagnostics")
    if not isinstance(diagnostic, dict) or diagnostic.get("passed") is not True:
        return issues + [f"{label} diagnostic evidence is incomplete or failed"]
    fields = (
        (("free_residual_relative", 1e-8), ("symmetry_error", 1e-12),
         ("energy_error_relative", 1e-8), ("resultant_balance_relative", 1e-7),
         ("moment_balance_relative", 1e-7)) if independent else
        (("residual_max", 1e-8), ("symmetry_error", 1e-12),
         ("energy_error_max", 1e-8), ("element_energy_error_max", 1e-8),
         ("force_balance_max", 1e-7), ("moment_balance_max", 1e-7)))
    for key, limit in fields:
        value = _maximum(diagnostic.get(key))
        if value is None or value < 0 or value > limit:
            issues.append(f"{label} {key} is missing, invalid, or exceeds {limit:g}")
    pivot_key = "lu_pivot_min_max_ratio" if independent else "minimum_relative_lu_pivot"
    pivot = diagnostic.get(pivot_key)
    if not _finite(pivot, positive=True) or pivot < 1e-14:
        issues.append(f"{label} LU pivot diagnostic is missing or indicates an unresolved near mechanism")
    if row.get("ill_conditioned") or diagnostic.get("ill_conditioned"):
        issues.append(f"{label} is explicitly flagged ill-conditioned")
    gram = row.get("load_gram")
    valid_gram = (isinstance(gram, (list, tuple)) and len(gram) == 2
                  and all(isinstance(line, (list, tuple)) and len(line) == 2 for line in gram)
                  and all(_finite(value) for line in gram for value in line))
    if not valid_gram:
        issues.append(f"{label} raw two-load Gram evidence is missing or invalid")
    else:
        scale = math.hypot(*(float(value) for line in gram for value in line))
        reciprocity = math.sqrt(2.)*abs(float(gram[0][1])-float(gram[1][0]))/max(scale, 1e-300)
        if reciprocity > 1e-8:
            issues.append(f"{label} raw Gram reciprocity exceeds 1e-8")
    return issues


def _geometry_status(records, audit):
    """Conservative geometry classification without changing connectivity."""
    binary = bool(records) and records[0].get("representation") == "binary"
    if any(r.get("status") in ("geometry_failed", "geometry_rejected") for r in records):
        return "geometry_rejected", "An archived/new geometry check rejected this frozen geometry"
    for record in records:
        old = record.get("geometry_status")
        if isinstance(old, dict) and old.get("valid") is False:
            return "geometry_rejected", "A frozen geometry fails edge connectivity or pad preservation"
        if record.get("pads_preserved") is False:
            return "geometry_rejected", "Required preserved pads are not established"
    if audit is None:
        if binary:
            return "geometry_unassessed", "Phase 1.1 local point-contact audit has not been supplied"
        return "not_applicable_density", "A relaxed density model is not an accepted binary geometry"
    if not isinstance(audit, dict):
        return "geometry_unassessed", "Geometry audit must be a mapping"
    if records:
        for key, record_key in (("candidate_hash", "hash"), ("hash", "hash"),
                                ("source_geometry_hash", "source_geometry_hash")):
            if key in audit and audit[key] != records[0].get(record_key):
                return "geometry_unassessed", "Geometry audit identity does not match the frozen records"
    # A contact audit may distinguish global edge connectivity from the local
    # diagonal-contact issue. Disconnected material is rejected before FEM.
    if (audit.get("edge_connected") is False or audit.get("pads_preserved") is False
            or audit.get("geometry_valid") is False
            or audit.get("geometry_rejected") is True
            or audit.get("status") in ("geometry_failed", "geometry_rejected")):
        return "geometry_rejected", "Independent geometry audit rejects connectivity or pads"
    contact_count = audit.get("local_point_contact_count", audit.get("vertex_only_contact_count", 0))
    if not _finite(contact_count) or contact_count < 0 or int(contact_count) != contact_count:
        return "geometry_unassessed", "Local point-contact count is invalid"
    if contact_count > 0 or audit.get("status") == "unresolved_point_contacts":
        return "unresolved_point_contacts", (
            "Globally connected material has local vertex-only contacts; shared displacement nodes "
            "impose a contact interpretation that numerical refinement alone does not validate")
    if audit.get("valid") is False:
        return "geometry_rejected", "Independent geometry audit explicitly rejects this geometry"
    checked_contacts = (audit.get("local_point_contact_detection_checked") is True
                        or audit.get("local_point_contacts_checked") is True)
    if (checked_contacts
            and (audit.get("valid") is True or audit.get("edge_connected") is True)
            and audit.get("pads_preserved") is True):
        return "geometry_valid", "Edge connectivity, preserved pads and absence of local point contacts established"
    return "geometry_unassessed", "Geometry audit lacks explicit connectivity, pad or local-contact evidence"


def _constraint_status(records):
    known = [r for r in records if "volume_fraction" in r or "volume_budget_fraction" in r]
    if not known:
        return "unknown", "Physical volume and budget are not recorded"
    if any(not _finite(r.get("volume_fraction")) or r["volume_fraction"] < 0
           or not _finite(r.get("volume_budget_fraction"), positive=True)
           or r["volume_budget_fraction"] > 1 for r in known):
        return "unknown", "Invalid physical-volume or budget evidence"
    budgets = [float(r["volume_budget_fraction"]) for r in known]
    volumes = [float(r["volume_fraction"]) for r in known]
    if max(budgets)-min(budgets) > 1e-12 or max(volumes)-min(volumes) > 1e-10:
        return "unknown", "Frozen volume or budget is inconsistent across records"
    # Full-volume normalized tolerance is exactly the original 1e-6 V_Omega.
    if any(volume > budget+1e-6 for volume, budget in zip(volumes, budgets)):
        return "over_budget", "Measured material volume exceeds the declared budget and tolerance"
    return "within_budget", "Measured material volume satisfies the declared budget and tolerance"


def reference_label(records, *, planned_levels=None, historical_label=None, geometry_audit=None):
    """Return a new label; preserve all input records and historical evidence.

    Rows use the existing ``evaluate_candidate`` schema. Failed/preflight rows
    may be incomplete: they produce an unresolved result, never an exception or
    an invented metric. Numerical acceptance and budget/geometry eligibility are
    separate. ``response_angle_degrees`` defaults to zero for historical rows;
    all selected and historical rows must declare the same angle. ``passed``
    means numerical gates for that angle only; use
    ``eligible_reference`` when selecting a geometrically and mass eligible
    binary reference. A positive numerical label remains an approximation to
    the prescribed continuum model, not exact physical ground truth.
    """
    result = {
        "schema_version": "phase1.1-reference-label-v1",
        "numerical_status": "unresolved", "passed": False,
        "eligible_reference": False, "classification": "unresolved",
        "verified_metrics": ["nominal_compliance_n_mm", "nominal_max_pad_displacement_mm"],
        "physical_validation": "not_assessed", "finest_is_exact_truth": False,
        "historical_label": deepcopy(historical_label), "historical_label_modified": False,
        "reason_codes": [], "reasons": [], "selected_levels_mm": [],
        "reference_estimate": None,
        "thresholds": {"compliance_relative": .01, "pad_relative": .02,
                       "compliance_monotonicity_relative": _MONOTONICITY_TOLERANCE,
                       "compliance_trend_comparison_slack": _TREND_TOLERANCE,
                       "volume_full_domain_fraction_slack": 1e-6,
                       "same_h_q4_p2_agreement_required": True},
    }
    try:
        records = list(records)
    except TypeError:
        records = []
        result["reason_codes"].append("invalid_records")
        result["reasons"].append("Records must be an iterable of mappings")
    malformed = any(not isinstance(r, dict) for r in records)
    if malformed:
        records = [r for r in records if isinstance(r, dict)]
        result["reason_codes"].append("invalid_records")
        result["reasons"].append("At least one record is not a mapping")
    geometry, geometry_reason = _geometry_status(records, geometry_audit)
    constraint, constraint_reason = _constraint_status(records)
    result.update(geometry_status=geometry, geometry_reason=geometry_reason,
                  constraint_status=constraint, constraint_reason=constraint_reason)

    def issue(code, detail):
        if code not in result["reason_codes"]:
            result["reason_codes"].append(code)
        result["reasons"].append(detail)

    def finish():
        if geometry == "geometry_rejected":
            result["classification"] = "geometry_rejected"
        elif constraint == "over_budget":
            result["classification"] = "constraint_over_budget"
        elif result["passed"] and geometry == "geometry_valid" and constraint == "within_budget":
            result["classification"] = "numerically_accepted"
            result["eligible_reference"] = True
        else:
            result["classification"] = "unresolved"
        result["status"] = result["classification"]
        return result

    if not records:
        issue("missing_records", "No mechanics records are available")
        return finish()
    for r in records:
        if any(not isinstance(r.get(k), str) or not r[k] for k in _STRING_IDENTITY):
            issue("incomplete_provenance", "Every frozen record must carry candidate, field, representation and problem/geometry identity")
        if not _finite(r.get("source_h_mm"), positive=True):
            issue("incomplete_provenance", "Every frozen record must carry finite positive source spacing")
        if not _finite(r.get("mesh_h_mm"), positive=True):
            issue("invalid_mesh_levels", "Every attempt, including a resource stop, needs finite positive mesh spacing")
    if result["reason_codes"]:
        return finish()
    identity = tuple(records[0][key] for key in _IDENTITY)
    if any(tuple(r[key] for key in _IDENTITY) != identity for r in records[1:]):
        issue("mixed_provenance", "Cannot combine different frozen candidates, geometries, representations or problems")
        return finish()
    result["frozen_identity"] = dict(zip(_IDENTITY, identity))
    angles = [r.get("response_angle_degrees", 0.) for r in records]
    if not all(_finite(angle) for angle in angles):
        issue("invalid_response_setting", "Response angles must be finite numbers in degrees")
        return finish()
    if any(float(angle) != float(angles[0]) for angle in angles[1:]):
        issue("mixed_response_setting", "Cannot combine response metrics from different load angles")
        return finish()
    angle = float(angles[0])
    result["frozen_identity"]["response_angle_degrees"] = angle
    result["response_angle_degrees"] = angle
    result["verified_metrics"] = ["compliance_n_mm", "max_pad_displacement_mm"]
    result["response_scope"] = f"Fixed {angle:g}-degree load; no other angle, angular maximum, stress or physical claim"
    source_h = float(records[0]["source_h_mm"])
    attempted = [float(r["mesh_h_mm"]) for r in records]
    result["attempted_levels_mm"] = sorted(attempted, reverse=True)
    if len(set(attempted)) != len(attempted):
        issue("duplicate_mesh_level", "Repeated mesh-level rows need a separately identified replay; do not choose a favorable duplicate")
        return finish()
    if planned_levels is None:
        levels = sorted(attempted, reverse=True)
    else:
        try:
            raw_levels = list(planned_levels)
        except TypeError:
            raw_levels = []
        if not raw_levels or not all(_finite(h, positive=True) for h in raw_levels):
            issue("invalid_planned_levels", "Planned levels must be a nonempty finite positive sequence")
            return finish()
        levels = sorted([float(h) for h in raw_levels], reverse=True)
        if len(set(levels)) != len(levels):
            issue("invalid_planned_levels", "Planned levels contain duplicates")
            return finish()
        # Caller may specify only the rolling window, but may never erase a
        # more recent finer attempt merely because that attempt failed.
        if any(h < min(levels) for h in attempted):
            issue("finer_attempt_omitted", "Planned levels omit a finer recorded attempt")
            return finish()
    result["committed_levels_mm"] = levels
    if len(levels) < 3:
        issue("fewer_than_three_levels", "Three nested refinement levels are required")
        return finish()
    selected = levels[-3:]
    result["selected_levels_mm"] = selected
    if any(not math.isclose(a/b, 2., rel_tol=0., abs_tol=1e-12) for a, b in zip(selected[:-1], selected[1:])):
        issue("nonnested_levels", "The final three spacings must be successive factors of two")
    if any(source_h/h < 1 or not math.isclose(source_h/h, round(source_h/h), rel_tol=0., abs_tol=1e-10) for h in selected):
        issue("source_grid_mismatch", "Each selected mesh must be an integer subdivision of the unchanged source grid")
    by_h = {float(r["mesh_h_mm"]): r for r in records}
    if any(h not in by_h for h in selected):
        issue("missing_committed_level", "A final committed refinement level has no result; older successful values cannot replace it")
    if result["reason_codes"]:
        return finish()
    selected_rows = [by_h[h] for h in selected]
    for r in selected_rows:
        h = r["mesh_h_mm"]
        if r.get("status") != "evaluated":
            issue("unsuccessful_selected_attempt", f"h={h:g} has status {r.get('status')!r}; no accepted solve is inferred")
        if r.get("pads_preserved") is not True:
            issue("pads_unverified", f"h={h:g} does not establish the complete preserved pads")
        independent = r.get("independent")
        if not isinstance(independent, dict):
            issue("independent_missing", f"h={h:g} has no independent response record")
            continue
        for label, part in (("Q4", r), ("P2", independent)):
            if any(not _finite(part.get(key), positive=True) for key in _METRICS):
                issue("invalid_metrics", f"h={h:g} {label} nominal metrics are missing, nonfinite or nonpositive")
            if label == "P2" and "mesh_h_mm" in part and part["mesh_h_mm"] != h:
                issue("independent_mesh_mismatch", f"h={h:g} has an independent result from another spacing")
            for detail in _diagnostic_issues(part, independent=label == "P2"):
                issue("solver_diagnostics_failed", f"h={h:g}: {detail}")
    if result["reason_codes"]:
        return finish()

    def change(a, b, key):
        return abs(float(a[key])-float(b[key]))/abs(float(b[key]))

    def trends(sequence):
        coarse, middle, fine = sequence
        previous = change(coarse, middle, "compliance_n_mm")
        last = change(middle, fine, "compliance_n_mm")
        pad = change(middle, fine, "max_pad_displacement_mm")
        values = [float(r["compliance_n_mm"]) for r in sequence]
        monotone = all(a <= b+_MONOTONICITY_TOLERANCE*max(abs(a), abs(b), 1e-12)
                       for a, b in zip(values[:-1], values[1:]))
        return {"previous_compliance_change": previous, "last_compliance_change": last,
                "last_pad_change": pad,
                "compliance_changes_decrease": last <= previous+_TREND_TOLERANCE,
                "compliance_nondecreasing": monotone}

    q4 = trends(selected_rows)
    p2 = trends([r["independent"] for r in selected_rows])
    fine = selected_rows[-1]
    agreement_c = change(fine, fine["independent"], "compliance_n_mm")
    agreement_pad = change(fine, fine["independent"], "max_pad_displacement_mm")
    gates = {
        "q4_compliance_refinement": q4["last_compliance_change"] < .01,
        "p2_compliance_refinement": p2["last_compliance_change"] < .01,
        "q4_pad_refinement": q4["last_pad_change"] < .02,
        "p2_pad_refinement": p2["last_pad_change"] < .02,
        "q4_decreasing_compliance_changes": q4["compliance_changes_decrease"],
        "p2_decreasing_compliance_changes": p2["compliance_changes_decrease"],
        "q4_compliance_monotonicity": q4["compliance_nondecreasing"],
        "p2_compliance_monotonicity": p2["compliance_nondecreasing"],
        "finest_same_h_compliance_agreement": agreement_c < .01,
        "finest_same_h_pad_agreement": agreement_pad < .02,
    }
    result.update(q4_refinement=q4, p2_refinement=p2, gates=gates,
                  finest_compliance_agreement=agreement_c,
                  finest_pad_agreement=agreement_pad,
                  finest_observed_metrics={"mesh_h_mm": selected[-1],
                                           "Q4": {k: float(fine[k]) for k in _METRICS},
                                           "P2": {k: float(fine["independent"][k]) for k in _METRICS}})
    if all(gates.values()):
        result["passed"] = True
        result["numerical_status"] = "numerically_accepted"
        result["reference_estimate"] = {
            "kind": "accepted_numerical_estimate_not_exact_truth",
            "mesh_h_mm": selected[-1], "element": "P2_triangle",
            **{k: float(fine["independent"][k]) for k in _METRICS},
        }
    else:
        issue("resolution_gates_failed", "Failed gates: "+", ".join(k for k, passed in gates.items() if not passed))
    return finish()
