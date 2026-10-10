"""Deterministic diagnostics for frozen pixel reconstructions, without FEA.

All geometry remains a union of origin-aligned square cells. Axis-run lengths
describe pixel spans; they are NOT minimum Euclidean ligaments, feature-size
guarantees, manufacturing constraints, or strength criteria. This module does
not smooth, repair, remove islands, refilter, or write any artifact.
"""
from collections import Counter, defaultdict
from dataclasses import asdict
import hashlib
import json
from time import perf_counter, process_time

import numpy as np

from .geometry import bracket_pads
from .models import Material


THRESHOLDS = (.35, .40, .45, .50, .55, .60, .65)
SCHEMA_VERSION = "prometheus.binary-audit.v1"


def _json_hash(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _regular_grid(mesh):
    """Validate the geometric contract used by exact cell counting and runs."""
    indices = np.asarray(mesh.cell_indices)
    if len(np.unique(indices, axis=0)) != len(indices):
        raise ValueError("Cell grid indices must be unique")
    expected = indices[:, None, :]*mesh.h + mesh.h*np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    if not np.allclose(mesh.nodes[mesh.cells], expected, rtol=0, atol=1e-10):
        raise ValueError("Binary audit requires origin-aligned regular square cells")
    return indices


def _binary(mesh, binary):
    values = np.asarray(binary)
    if values.shape != (len(mesh.cells),) or not np.isin(values, [0, 1]).all():
        raise ValueError("Supply one binary zero/one value per cell")
    return values.astype(bool, copy=True)


def _density(mesh, physical):
    values = np.asarray(physical, dtype=float)
    if (values.shape != (len(mesh.cells),) or not np.isfinite(values).all()
            or np.any(values < 0) or np.any(values > 1)):
        raise ValueError("Density must be a finite [0,1] value per cell")
    return values


def _mask(mesh, values, label):
    result = np.asarray(values)
    if result.dtype != bool or result.shape != (len(mesh.cells),):
        raise ValueError(label+" must be a Boolean value per cell")
    return result.copy()


def binary_geometry_hash(mesh, binary):
    """Hash sorted retained square coordinates, independent of array ordering.

    The hash is for the declared pixel tessellation: equivalent refinement into
    smaller squares is not canonicalized. The Phase 1.1 catalogue freezes h=2.
    """
    indices = _regular_grid(mesh)
    solid = _binary(mesh, binary)
    selected = indices[solid]
    selected = selected[np.lexsort((selected[:, 0], selected[:, 1]))]
    return _json_hash({"representation": "union_of_grid_squares",
                       "cell_size_mm": float(mesh.h), "solid_cell_indices": selected.tolist()})


def source_fingerprint(mesh, physical, source_input_hash):
    """Record caller-supplied archive-byte digest and independently hash values.

    ``source_input_hash`` must be the SHA-256 of the archived input bytes (or a
    manifest explicitly identifying them). This pure helper cannot authenticate
    the supplied digest against a file; archive loading is the caller's job.
    """
    _regular_grid(mesh)
    values = _density(mesh, physical)
    if (not isinstance(source_input_hash, str) or len(source_input_hash) != 64
            or any(c not in "0123456789abcdef" for c in source_input_hash)):
        raise ValueError("source_input_hash must be a lowercase SHA-256 hex digest")
    mesh_hash = _json_hash({"nodes": mesh.nodes.tolist(), "cells": mesh.cells.tolist(),
                           "cell_indices": mesh.cell_indices.tolist(), "h_mm": float(mesh.h)})
    return {"source_input_sha256": source_input_hash,
            "source_input_hash_origin": "caller_supplied_archive_digest",
            "source_field_sha256": hashlib.sha256(values.astype("<f8").tobytes()).hexdigest(),
            "source_mesh_sha256": mesh_hash, "source_h_mm": float(mesh.h)}


def _components(indices, solid):
    remaining = {tuple(map(int, index)) for index in indices[solid]}
    components = []
    while remaining:
        start = min(remaining, key=lambda cell: (cell[1], cell[0]))
        remaining.remove(start)
        pending, component = [start], []
        while pending:
            x, y = pending.pop()
            component.append((x, y))
            for neighbor in ((x-1, y), (x+1, y), (x, y-1), (x, y+1)):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    pending.append(neighbor)
        components.append(sorted(component, key=lambda cell: (cell[1], cell[0])))
    return components


def _axis_runs(indices, solid, h, axis):
    """Each maximal contiguous run counts once, not once per constituent cell."""
    lines = defaultdict(list)
    for index in indices[solid]:
        lines[int(index[1-axis])].append(int(index[axis]))
    lengths = []
    for key in sorted(lines):
        along = sorted(lines[key])
        length = 1
        for previous, current in zip(along[:-1], along[1:]):
            if current == previous+1:
                length += 1
            else:
                lengths.append(length)
                length = 1
        lengths.append(length)
    counts = Counter(lengths)
    physical = np.asarray(lengths, dtype=float)*h
    return {"run_count": len(lengths),
            "minimum_mm": float(physical.min()) if lengths else None,
            "median_mm": float(np.median(physical)) if lengths else None,
            "maximum_mm": float(physical.max()) if lengths else None,
            "histogram": [{"length_cells": int(length), "length_mm": float(length*h),
                           "run_count": int(count)} for length, count in sorted(counts.items())],
            "weighting": "one count per maximal contiguous axis-aligned run"}


def _pixel_corners(indices, solid, h):
    retained = {tuple(map(int, index)) for index in indices[solid]}
    vertices = set()
    for x, y in retained:
        vertices.update(((x, y), (x+1, y), (x+1, y+1), (x, y+1)))
    convex, reentrant, vertex_only = [], [], []
    for x, y in sorted(vertices, key=lambda p: (p[1], p[0])):
        occupied = [(x-1, y-1) in retained, (x, y-1) in retained,
                    (x, y) in retained, (x-1, y) in retained]
        count = sum(occupied)
        point = [float(x*h), float(y*h)]
        if count == 1:
            convex.append(point)
        elif count == 3:
            reentrant.append(point)
        elif count == 2 and occupied[0] == occupied[2] and occupied[1] == occupied[3]:
            vertex_only.append(point)
    return {"convex_90_degree_count": len(convex), "reentrant_270_degree_count": len(reentrant),
            "vertex_only_contact_count": len(vertex_only), "convex_vertices_mm": convex,
            "reentrant_vertices_mm": reentrant, "vertex_only_contacts_mm": vertex_only,
            "interpretation": "pixel boundary geometry only; no stress or fillet inference"}


def binary_diagnostics(mesh, binary, material=Material(), budget_fraction=None,
                       reference_density=None, support_mask=None, load_mask=None):
    """Return JSON-ready geometry, exact cell measures, and bounded diagnostics.

    Default support/load masks are the canonical bracket pads, whose required
    areas are 400/200 mm². Explicit masks permit analytical rectangular fixtures.
    ``offending_area_mm2`` is material outside the first support-pad component,
    matching the Phase 1 connectivity convention when the full pads are present.
    """
    indices = _regular_grid(mesh)
    solid = _binary(mesh, binary)
    # Revalidate user-supplied material records even when using duck typing.
    if (not np.isfinite([material.thickness, material.density]).all()
            or min(material.thickness, material.density) <= 0):
        raise ValueError("Thickness and mass density must be positive and finite")
    if budget_fraction is not None and (not np.isfinite(budget_fraction) or not 0 < budget_fraction <= 1):
        raise ValueError("Budget fraction must lie in (0,1]")
    centers = mesh.centroids
    pads = bracket_pads(mesh)
    support = (pads & (centers[:, 1] >= 70)) if support_mask is None else _mask(mesh, support_mask, "Support mask")
    load = (pads & (centers[:, 0] >= 110)) if load_mask is None else _mask(mesh, load_mask, "Load mask")
    area_per_cell = float(mesh.h**2)
    expected_support_area = 400. if support_mask is None else area_per_cell*int(support.sum())
    expected_load_area = 200. if load_mask is None else area_per_cell*int(load.sum())
    support_present = bool(support.any() and np.all(solid[support]) and
                           np.isclose(support.sum()*area_per_cell, expected_support_area, rtol=0, atol=1e-8))
    load_present = bool(load.any() and np.all(solid[load]) and
                        np.isclose(load.sum()*area_per_cell, expected_load_area, rtol=0, atol=1e-8))
    support_cells = {tuple(map(int, index)) for index in indices[support]}
    load_cells = {tuple(map(int, index)) for index in indices[load]}
    components = _components(indices, solid)
    rows = []
    for component_id, component in enumerate(components):
        cells = set(component)
        points = np.asarray(component)
        rows.append({"component_id": component_id, "cell_count": len(component),
                     "area_mm2": len(component)*area_per_cell,
                     "support_pad_cell_count": len(cells & support_cells),
                     "load_pad_cell_count": len(cells & load_cells),
                     "bounds_mm": [float(points[:, 0].min()*mesh.h), float(points[:, 1].min()*mesh.h),
                                   float((points[:, 0].max()+1)*mesh.h), float((points[:, 1].max()+1)*mesh.h)]})
    anchors = [r["component_id"] for r in rows if r["support_pad_cell_count"]]
    anchor = anchors[0] if anchors else None
    offenders = [cell for i, component in enumerate(components) if i != anchor for cell in component]
    pads_connected = bool(support_present and load_present and anchor is not None
                          and support_cells | load_cells <= set(components[anchor]))
    connected = len(components) == 1
    area = int(solid.sum())*area_per_cell
    total_area = len(solid)*area_per_cell
    volume = area*material.thickness
    reference_area = None if reference_density is None else float(_density(mesh, reference_density).sum()*area_per_cell)
    budget_area = None if budget_fraction is None else float(budget_fraction*total_area)
    delta = None if budget_area is None else area-budget_area
    density_delta = None if reference_area is None else area-reference_area
    x_runs = _axis_runs(indices, solid, mesh.h, 0)
    y_runs = _axis_runs(indices, solid, mesh.h, 1)
    return {"schema_version": SCHEMA_VERSION, "geometry_sha256": binary_geometry_hash(mesh, solid),
            "cell_size_mm": float(mesh.h), "domain_area_mm2": total_area,
            "solid_cell_count": int(solid.sum()), "solid_area_mm2": area,
            "material_volume_mm3": volume, "mass_kg": volume*material.density,
            "volume_fraction": area/total_area, "budget_fraction": budget_fraction,
            "budget_area_mm2": budget_area, "area_minus_budget_mm2": delta,
            "volume_minus_budget_mm3": None if delta is None else delta*material.thickness,
            "mass_minus_budget_kg": None if delta is None else delta*material.thickness*material.density,
            "budget_feasible": None if delta is None else bool(delta <= 1e-6*total_area),
            "budget_area_tolerance_mm2": 1e-6*total_area,
            "reference_density_area_mm2": reference_area,
            "area_change_from_density_mm2": density_delta,
            "mass_change_from_density_kg": None if density_delta is None else density_delta*material.thickness*material.density,
            "area_method": "retained square count times h squared; no contour approximation",
            "edge_component_count": len(components), "edge_connected": connected,
            "components": rows, "support_pad_preserved": support_present,
            "load_pad_preserved": load_present, "pads_connected": pads_connected,
            "all_retained_material_support_connected": bool(connected and pads_connected),
            "geometry_valid": bool(connected and pads_connected),
            "support_component_id": anchor, "offending_area_mm2": len(offenders)*area_per_cell,
            "offending_cell_indices": [list(cell) for cell in offenders],
            "axis_runs": {"x": x_runs, "y": y_runs},
            "minimum_axis_run_mm": min(x_runs["minimum_mm"], y_runs["minimum_mm"]) if area else None,
            "minimum_ligament_mm": None, "minimum_ligament_status": "not_measured",
            "axis_run_interpretation": "grid-aligned maximal spans, not Euclidean ligament thickness or minimum-feature guarantees",
            "pixel_corners": _pixel_corners(indices, solid, mesh.h),
            "physical_validation": "not_assessed", "manufacturability": "not_assessed",
            "factorizations": 0, "rhs_solves": 0}


def threshold_reconstruction(mesh, physical, threshold, source_input_hash, material=Material(),
                             budget_fraction=None, aliases=None):
    """Generate a new threshold candidate, preserving every canonical pad cell.

    Candidate identity includes source/input provenance, rule/threshold, and
    output geometry. Geometry identity allows equivalent outputs to reuse an
    existing analysis identity through an explicit alias registry. Registry and
    input arrays are never modified. Returned binary arrays are read-only.
    """
    started, cpu = perf_counter(), process_time()
    values = _density(mesh, physical)
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("Threshold must be finite and lie in [0,1]")
    provenance = source_fingerprint(mesh, values, source_input_hash)
    pads = bracket_pads(mesh)
    solid = values >= threshold
    overridden = int(np.count_nonzero(pads & ~solid))
    solid[pads] = True
    diagnostics = binary_diagnostics(mesh, solid, material, budget_fraction, values)
    rule = {"name": "physical_density_ge_threshold_then_solid_pad_override",
            "threshold": float(threshold), "threshold_hex": float(threshold).hex(),
            "pad_mask_sha256": hashlib.sha256(pads.astype(np.uint8).tobytes()).hexdigest(),
            "repair": "none", "island_removal": "none"}
    geometry_hash = diagnostics["geometry_sha256"]
    identity = {"schema_version": SCHEMA_VERSION, **provenance, "rule": rule,
                "geometry_sha256": geometry_hash}
    candidate_id = "binary-"+_json_hash(identity)
    alias = None if aliases is None else aliases.get(geometry_hash)
    if alias is not None and (not isinstance(alias, str) or not alias):
        raise ValueError("Geometry aliases must identify an existing candidate by a nonempty string")
    solid.setflags(write=False)
    record = {**identity, "candidate_id": candidate_id,
              "geometry_id": "pixel-geometry-"+geometry_hash,
              "geometry_alias_of": alias if alias != candidate_id else None,
              "material": asdict(material), "pad_cells_overridden": overridden,
              "diagnostics": diagnostics, "generation_wall_seconds": perf_counter()-started,
              "generation_cpu_seconds": process_time()-cpu,
              "factorizations": 0, "rhs_solves": 0, "status": "generated_not_analyzed"}
    return solid, record


def threshold_sweep(mesh, physical, source_input_hash, material=Material(), budget_fraction=None,
                    thresholds=THRESHOLDS, aliases=None):
    """Generate the predefined threshold catalogue in the declared order."""
    thresholds = tuple(thresholds)
    if not thresholds or not np.isfinite(thresholds).all() or len(set(thresholds)) != len(thresholds):
        raise ValueError("Threshold sweep needs unique finite thresholds")
    registry = {} if aliases is None else dict(aliases)
    generated = []
    for threshold in thresholds:
        binary, record = threshold_reconstruction(mesh, physical, threshold, source_input_hash,
                                                  material, budget_fraction, registry)
        registry.setdefault(record["geometry_sha256"], record["candidate_id"])
        generated.append((binary, record))
    return generated
