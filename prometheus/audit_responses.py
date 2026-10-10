"""Additive Phase 1.1 response projection and topology-based preflight.

Projection reuses already solved right/down basis fields. It neither solves a
new equilibrium system nor upgrades an acceptance label. Resource predictions
count the unchanged shared-node pixel mesh; they do not estimate factor fill,
prove a memory bound, or imply that a matrix was assembled.
"""
from copy import deepcopy
import hashlib
import json
import math

import numpy as np


_DIRECTIONAL_METRICS = (
    "compliance_n_mm", "max_pad_displacement_mm", "load_weighted_displacement_mm",
    "max_displacement_all_nodes_mm", "max_material_displacement_mm",
    "small_deformation_review_trigger",
)
_MAX_P2_DOFS = 1_500_000
_MAX_P2_COO_BYTES = 4*1024**3


def _project_part(part, arrays, node_key, displacement_key, coefficients):
    nodes = np.asarray(arrays[node_key], dtype=float)
    basis = np.asarray(arrays[displacement_key], dtype=float)
    gram = np.asarray(part["load_gram"], dtype=float)
    if (nodes.ndim != 2 or nodes.shape[1] != 2 or not len(nodes)
            or not np.isfinite(nodes).all() or basis.shape != (2*len(nodes), 2)
            or not np.isfinite(basis).all() or gram.shape != (2, 2)
            or not np.isfinite(gram).all()):
        raise ValueError("Response projection requires finite nodes, (ndof,2) basis displacements and a 2x2 load Gram matrix")
    x, y = nodes.T
    tol = 1e-10
    pads = (((x >= -tol) & (x <= 40+tol) & (y >= 70-tol) & (y <= 80+tol))
            | ((x >= 110-tol) & (x <= 120+tol) & (y >= 10-tol) & (y <= 30+tol)))
    if not pads.any():
        raise ValueError("No nodes lie on the declared canonical pads")
    displacement = (basis@coefficients).reshape(-1, 2)
    magnitude = np.linalg.norm(displacement, axis=1)
    return {
        "compliance_n_mm": float(coefficients@gram@coefficients),
        "max_pad_displacement_mm": float(magnitude[pads].max()),
        "max_displacement_all_nodes_mm": float(magnitude.max()),
        "small_deformation_review_trigger": bool(magnitude[pads].max()/120 > .01),
    }


def project_response(record, arrays, angle_degrees):
    """Deep-copy a mechanics record and calculate its declared-angle response.

    The physical bases are RIGHT and DOWN, including their original resultant
    magnitude. Nominal records need no signature change; nonzero angles have a
    stable signature derived from the unchanged base signature and angle.
    Raw basis Gram matrices and solver diagnostics are retained. Directional
    original values are preserved in ``base_response_metrics`` rather than
    being mislabeled as responses at another angle when a projection cannot be
    made. No missing or unsuccessful solve is supplied with invented metrics.
    """
    if not isinstance(record, dict):
        raise ValueError("A response record must be a mapping")
    if isinstance(angle_degrees, (bool, np.bool_)):
        raise ValueError("Response angle must be a finite number in degrees")
    try:
        angle = float(angle_degrees)
    except (ValueError, TypeError) as exc:
        raise ValueError("Response angle must be a finite number in degrees") from exc
    if not math.isfinite(angle):
        raise ValueError("Response angle must be finite")
    if angle == 0:
        angle = 0.  # normalize negative zero for stable identity
    result = deepcopy(record)
    base = record.get("base_problem_signature", record.get("problem_signature"))
    result["base_problem_signature"] = base
    result["response_angle_degrees"] = angle
    if isinstance(base, str) and base:
        signature = base if angle == 0 else hashlib.sha256(json.dumps(
            {"base_problem_signature": base, "response_angle_degrees": angle},
            sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
        result["problem_signature"] = signature
    result["response_projection_factorizations"] = 0
    result["response_projection_rhs_solves"] = 0
    result["response_projection_status"] = "not_projected"
    # Projection cannot transfer the original nominal verification flag.
    result["base_independent_verified"] = record.get("base_independent_verified", record.get("independent_verified", False))
    result["independent_verified"] = False
    for part in (result, result.get("independent")):
        if not isinstance(part, dict):
            continue
        if "base_response_metrics" not in part:
            part["base_response_metrics"] = {key: deepcopy(part[key]) for key in _DIRECTIONAL_METRICS if key in part}
        for key in _DIRECTIONAL_METRICS:
            part.pop(key, None)
    if record.get("status") != "evaluated":
        result["response_projection_status"] = "not_projected_unsuccessful_attempt"
        return result
    if not isinstance(base, str) or not base:
        result["status"] = "response_projection_failed"
        result["response_projection_status"] = "missing_base_problem_signature"
        return result
    if not isinstance(arrays, dict):
        # NpzFile implements a Mapping-like files/keys interface without being
        # a dict; accept it without forcing a copy of every large field array.
        if not hasattr(arrays, "__contains__") or not hasattr(arrays, "__getitem__"):
            result["status"] = "response_projection_failed"
            result["response_projection_status"] = "missing_basis_fields"
            return result
    required = ("nodes", "u_basis")
    if any(key not in arrays for key in required) or "load_gram" not in record:
        result["status"] = "response_projection_failed"
        result["response_projection_status"] = "missing_q4_basis_fields"
        return result
    coefficients = np.array([math.sin(math.radians(angle)), math.cos(math.radians(angle))])
    result.update(_project_part(record, arrays, "nodes", "u_basis", coefficients))
    # Material-region displacement is optional, but never inherit a nominal
    # maximum under a new angle. Reconstruct it when its original data exist.
    if "cells" in arrays and "rho" in arrays:
        cells = np.asarray(arrays["cells"])
        rho = np.asarray(arrays["rho"], dtype=float)
        if (cells.ndim != 2 or cells.shape[1] != 4 or rho.shape != (len(cells),)
                or not np.issubdtype(cells.dtype, np.integer) or not np.isfinite(rho).all()
                or np.any(cells < 0) or np.any(cells >= len(arrays["nodes"]))):
            raise ValueError("Material displacement requires valid cells and per-cell densities")
        material_nodes = np.unique(cells[rho >= .9])
        magnitude = np.linalg.norm((np.asarray(arrays["u_basis"])@coefficients).reshape(-1, 2), axis=1)
        result["max_material_displacement_mm"] = float(magnitude[material_nodes].max()) if len(material_nodes) else None
    # The common evaluator's force arrays permit recovering the resultant
    # magnitude, avoiding an undocumented hardcoded 1000 N conversion.
    if "forces_basis" in arrays:
        forces = np.asarray(arrays["forces_basis"], dtype=float)
        if forces.shape != np.asarray(arrays["u_basis"]).shape or not np.isfinite(forces).all():
            raise ValueError("Forces must match the two displacement basis columns")
        resultant = (forces@coefficients).reshape(-1, 2).sum(axis=0)
        load = float(np.linalg.norm(resultant))
        if load > 0:
            result["load_weighted_displacement_mm"] = result["compliance_n_mm"]/load
    independent = record.get("independent")
    if isinstance(independent, dict):
        if (independent.get("solver_passed") is not True or "load_gram" not in independent
                or any(key not in arrays for key in ("p2_nodes", "p2_u_basis"))):
            result["status"] = "response_projection_failed"
            result["response_projection_status"] = "independent_response_unavailable"
            result["independent"]["solver_passed"] = False
            return result
        result["independent"].update(_project_part(independent, arrays, "p2_nodes", "p2_u_basis", coefficients))
        result["independent"]["response_angle_degrees"] = angle
        # These display caches refer to the projected response. Acceptance
        # still recomputes the gaps directly in reference_labels.py.
        for key, cache in (("compliance_n_mm", "q4_relative_compliance_difference"),
                           ("max_pad_displacement_mm", "q4_relative_pad_difference")):
            scale = abs(result["independent"][key])
            result["independent"][cache] = abs(result[key]-result["independent"][key])/scale if scale else None
    result["response_projection_status"] = "projected_both_solvers" if isinstance(independent, dict) else "projected_q4_only"
    return result


def resource_estimate(source_mesh, binary, target_h):
    """Predict uniform-refinement topology and apply declared preflight caps.

    With k subdivisions per source edge: V'=V+(k-1)E+(k-1)^2N,
    E'=kE+2k(k-1)N, N'=k²N. Splitting each quad creates one diagonal
    midpoint, so P2 nodes = V'+E'+N'. Counts assume unchanged shared-node
    connectivity and do not certify domain/contact validity.
    """
    if isinstance(target_h, (bool, np.bool_)):
        raise ValueError("Target spacing must be finite and positive")
    try:
        target_h, source_h = float(target_h), float(source_mesh.h)
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError("Source and target spacings must be finite and positive") from exc
    if not math.isfinite(target_h) or not math.isfinite(source_h) or min(target_h, source_h) <= 0:
        raise ValueError("Source and target spacings must be finite and positive")
    ratio = source_h/target_h
    if ratio < 1 or not math.isclose(ratio, round(ratio), rel_tol=0., abs_tol=1e-10):
        raise ValueError("Target mesh must be an integer subdivision of the source")
    k = int(round(ratio))
    cells = np.asarray(source_mesh.cells)
    nodes = np.asarray(source_mesh.nodes, dtype=float)
    mask = np.asarray(binary)
    if (cells.ndim != 2 or cells.shape[1] != 4 or not np.issubdtype(cells.dtype, np.integer)
            or nodes.ndim != 2 or nodes.shape[1] != 2 or not np.isfinite(nodes).all()
            or np.any(cells < 0) or np.any(cells >= len(nodes))
            or mask.shape != (len(cells),) or not np.isin(mask, [0, 1]).all()):
        raise ValueError("Expected a valid quad mesh and one binary value per source cell")
    active = cells[mask.astype(bool)]
    n = len(active)
    v = len(np.unique(active))
    edges = active[:, [[0, 1], [1, 2], [2, 3], [3, 0]]].reshape(-1, 2)
    e = len(np.unique(np.sort(edges, axis=1), axis=0))
    refined_v = v+(k-1)*e+(k-1)**2*n
    refined_e = k*e+2*k*(k-1)*n
    refined_n = k*k*n
    p2_dofs = 2*(refined_v+refined_e+refined_n)
    estimated_bytes = 2*refined_n*144*40+p2_dofs*64
    reasons = []
    if not n:
        reasons.append("Empty material geometry cannot be analyzed")
    if p2_dofs > _MAX_P2_DOFS:
        reasons.append("Predicted P2 DOFs exceed the declared 1500000 ceiling")
    if estimated_bytes > _MAX_P2_COO_BYTES:
        reasons.append("Conservative P2 assembly allocation estimate exceeds 4 GiB")
    return {
        "source_h_mm": source_h, "target_h_mm": target_h, "subdivision_factor": k,
        "source_material_quads": n, "source_active_vertices": v, "source_active_quad_edges": e,
        "q4_cells": refined_n, "q4_vertices": refined_v, "q4_edges": refined_e,
        "p2_triangles": 2*refined_n, "p2_dofs": p2_dofs,
        "estimated_p2_coo_bytes": estimated_bytes,
        "max_p2_dofs": _MAX_P2_DOFS, "max_p2_coo_bytes": _MAX_P2_COO_BYTES,
        "allowed": not reasons, "reasons": reasons,
        "estimate_scope": "Predicted unchanged shared-node topology and conservative assembly allocation; no factorization fill or measured memory",
        "matrix_assembled": False, "factorization_memory_estimate_bytes": None,
    }
