"""Analytical pixel-geometry fixtures; no mechanics solve or archive mutation."""
import hashlib
import json

import numpy as np
import pytest

from prometheus.binary_audit import (THRESHOLDS, binary_diagnostics, binary_geometry_hash,
    source_fingerprint, threshold_reconstruction, threshold_sweep, retain_pad_component)
from prometheus.geometry import Mesh, make_rect_mesh, make_l_mesh, bracket_pads
from prometheus.models import Material
from prometheus.optimization import simple_l_design


INPUT_HASH = hashlib.sha256(b"declared analytical input fixture").hexdigest()


def test_rectangle_exact_measures_and_axis_run_distribution():
    mesh = make_rect_mesh(8, 4, 2)
    solid = np.ones(len(mesh.cells), dtype=bool)
    support = mesh.cell_indices[:, 0] == 0
    load = mesh.cell_indices[:, 0] == 3
    result = binary_diagnostics(mesh, solid, budget_fraction=.5,
                                reference_density=np.full(8, .25), support_mask=support, load_mask=load)
    assert result["solid_area_mm2"] == 32
    assert result["material_volume_mm3"] == 160
    assert result["mass_kg"] == pytest.approx(160*2.7e-6)
    assert result["area_minus_budget_mm2"] == 16
    assert result["volume_minus_budget_mm3"] == 80
    assert result["mass_minus_budget_kg"] == pytest.approx(80*2.7e-6)
    assert result["area_change_from_density_mm2"] == 24
    assert not result["budget_feasible"]
    assert result["geometry_valid"] and result["edge_component_count"] == 1
    assert result["axis_runs"]["x"]["histogram"] == [{"length_cells": 4, "length_mm": 8., "run_count": 2}]
    assert result["axis_runs"]["y"]["histogram"] == [{"length_cells": 2, "length_mm": 4., "run_count": 4}]
    assert result["minimum_axis_run_mm"] == 4
    assert result["minimum_ligament_mm"] is None
    assert result["minimum_ligament_status"] == "not_measured"
    assert result["pixel_corners"]["convex_90_degree_count"] == 4
    assert result["pixel_corners"]["reentrant_270_degree_count"] == 0
    assert result["local_point_contact_count"] == 0
    assert result["factorizations"] == result["rhs_solves"] == 0
    json.dumps(result, allow_nan=False)


def test_disconnected_diagonal_cells_and_contact_coordinates():
    mesh = make_rect_mesh(2, 2, 1)
    result = binary_diagnostics(mesh, [1, 0, 0, 1],
                                support_mask=np.array([True, False, False, False]),
                                load_mask=np.array([False, False, False, True]))
    assert result["edge_component_count"] == 2
    assert not result["edge_connected"] and not result["geometry_valid"]
    assert result["offending_area_mm2"] == 1
    assert result["local_point_contact_count"] == 1
    assert result["local_point_contacts"] == [{"vertex_mm": [1., 1.],
                                                "solid_cell_indices": [[0, 0], [1, 1]],
                                                "solid_cell_ids": [0, 3], "shared_area_mm2": 0.}]


def test_global_edge_connectivity_does_not_exclude_a_local_point_contact():
    mesh = make_rect_mesh(4, 4, 1)
    x, y = mesh.cell_indices.T
    # Outer loop joins both diagonal cells by an alternate edge-connected path.
    solid = (x == 0) | (x == 3) | (y == 0) | (y == 3) | ((x == 1) & (y == 1)) | ((x == 2) & (y == 2))
    result = binary_diagnostics(mesh, solid, support_mask=x == 0, load_mask=x == 3)
    assert result["edge_connected"] and result["geometry_valid"]
    assert result["local_point_contact_detection_checked"]
    assert result["local_point_contact_count"] == 1
    assert result["local_point_contacts"][0]["vertex_mm"] == [2., 2.]
    assert result["local_point_contacts"][0]["solid_cell_ids"] == [5, 10]


def test_l_pixel_has_one_reentrant_corner_and_no_point_contact():
    mesh = make_rect_mesh(2, 2, 1)
    result = binary_diagnostics(mesh, [1, 1, 1, 0])
    assert result["edge_component_count"] == 1
    assert result["pixel_corners"]["reentrant_vertices_mm"] == [[1., 1.]]
    assert result["pixel_corners"]["convex_90_degree_count"] == 5
    assert result["local_point_contact_count"] == 0


def test_empty_geometry_has_no_fabricated_zero_ligament():
    mesh = make_rect_mesh(2, 2, 1)
    result = binary_diagnostics(mesh, np.zeros(4))
    assert result["solid_area_mm2"] == 0 and result["edge_component_count"] == 0
    assert not result["geometry_valid"]
    assert result["minimum_axis_run_mm"] is None
    assert result["axis_runs"]["x"]["run_count"] == 0
    assert result["local_point_contact_count"] == 0


def test_missing_domain_cells_act_as_void_in_local_contact_check():
    original = make_rect_mesh(2, 2, 1)
    selection = np.array([0, 3])
    mesh = Mesh(original.nodes, original.cells[selection], original.cell_indices[selection], 1)
    result = binary_diagnostics(mesh, [1, 1])
    assert result["local_point_contact_count"] == 1
    assert result["local_point_contacts"][0]["solid_cell_ids"] == [0, 1]


def test_canonical_full_geometry_and_pad_regions_have_exact_area():
    mesh = make_l_mesh(2)
    result = binary_diagnostics(mesh, np.ones(len(mesh.cells)), budget_fraction=.4)
    assert result["solid_area_mm2"] == 6400
    assert result["material_volume_mm3"] == 32000
    assert result["mass_kg"] == pytest.approx(.0864)
    assert result["pads_preserved"] and result["pads_connected"]
    assert result["offending_area_mm2"] == 0
    assert result["geometry_valid"] and result["local_point_contact_count"] == 0
    assert result["pixel_corners"]["reentrant_vertices_mm"] == [[40., 40.]]


def test_threshold_sweep_mass_is_nonincreasing_and_pads_stay_solid():
    mesh = make_l_mesh(2)
    physical = np.linspace(0, 1, len(mesh.cells))
    original = physical.copy()
    generated = threshold_sweep(mesh, physical, INPUT_HASH, budget_fraction=.4)
    assert len(generated) == len(THRESHOLDS) == 7
    previous = np.ones(len(mesh.cells), dtype=bool)
    for threshold, (binary, record) in zip(THRESHOLDS, generated):
        assert np.all(binary <= previous)
        assert np.all(binary[bracket_pads(mesh)])
        assert record["rule"]["threshold"] == threshold
        assert record["rule"]["repair"] == "none"
        assert record["diagnostics"]["solid_area_mm2"] == 4*binary.sum()
        assert record["status"] == "generated_not_analyzed"
        assert record["factorizations"] == record["rhs_solves"] == 0
        assert record["generation_wall_seconds"] >= 0 and record["generation_cpu_seconds"] >= 0
        json.dumps(record, allow_nan=False)
        with pytest.raises(ValueError):
            binary[0] = False
        previous = binary
    np.testing.assert_array_equal(physical, original)


def test_threshold_uses_greater_equal_and_not_a_hidden_mass_match():
    mesh = make_rect_mesh(2, 2, 1)
    binary, record = threshold_reconstruction(mesh, [0, .5, .5, 1], .5, INPUT_HASH,
                                              budget_fraction=.5)
    np.testing.assert_array_equal(binary, [False, True, True, True])
    assert record["diagnostics"]["volume_fraction"] == .75
    assert not record["diagnostics"]["budget_feasible"]


def test_identical_inputs_reproduce_ids_but_changed_rule_or_source_gets_new_id():
    mesh = make_rect_mesh(2, 2, 1)
    values = np.array([.1, .4, .6, .9])
    _, first = threshold_reconstruction(mesh, values, .5, INPUT_HASH)
    _, repeated = threshold_reconstruction(mesh, values, .5, INPUT_HASH)
    _, other_rule = threshold_reconstruction(mesh, values, .55, INPUT_HASH)
    _, other_source = threshold_reconstruction(mesh, values, .5, "0"*64)
    _, other_geometry = threshold_reconstruction(mesh, values, .65, INPUT_HASH)
    assert first["candidate_id"] == repeated["candidate_id"]
    assert first["candidate_id"] != other_rule["candidate_id"]
    assert first["geometry_id"] == other_rule["geometry_id"]
    assert first["candidate_id"] != other_source["candidate_id"]
    assert first["geometry_id"] != other_geometry["geometry_id"]


def test_same_geometry_aliases_an_existing_identity_without_mutating_registry():
    mesh = make_rect_mesh(2, 2, 1)
    values = np.ones(4)
    geometry_hash = binary_geometry_hash(mesh, values)
    aliases = {geometry_hash: "archived-reference-identity"}
    before = aliases.copy()
    generated = threshold_sweep(mesh, values, INPUT_HASH, aliases=aliases)
    assert all(record["geometry_alias_of"] == "archived-reference-identity" for _, record in generated)
    assert len({record["candidate_id"] for _, record in generated}) == 7
    assert len({record["geometry_id"] for _, record in generated}) == 1
    assert aliases == before
    generated = threshold_sweep(mesh, values, INPUT_HASH)
    assert generated[0][1]["geometry_alias_of"] is None
    assert generated[-1][1]["geometry_alias_of"] == generated[0][1]["candidate_id"]


def test_geometry_identity_is_order_independent_and_distinguishes_equal_area_shapes():
    mesh = make_rect_mesh(2, 2, 1)
    binary = np.array([1, 0, 1, 0])
    order = np.array([3, 1, 0, 2])
    reordered = Mesh(mesh.nodes, mesh.cells[order], mesh.cell_indices[order], 1)
    assert binary_geometry_hash(mesh, binary) == binary_geometry_hash(reordered, binary[order])
    assert binary_geometry_hash(mesh, binary) != binary_geometry_hash(mesh, 1-binary)


def test_declared_island_removal_creates_new_geometry_provenance_and_exact_new_mass():
    mesh = make_l_mesh(2)
    connected = simple_l_design(mesh, .35)
    assert not connected[0]
    original = connected.copy()
    original[0] = True  # A distant lower-left island, not a pad or connecting path.
    before = binary_diagnostics(mesh, original)
    assert before["edge_component_count"] == 2 and before["pads_connected"]
    copied = original.copy()
    repaired, record = retain_pad_component(mesh, original, INPUT_HASH, Material(), .35)
    np.testing.assert_array_equal(original, copied)
    np.testing.assert_array_equal(repaired, connected)
    assert record["rule"]["name"] == "retain_pad_connected_component"
    assert record["rule"]["removed_cell_indices"] == [[0, 0]]
    assert record["geometry_changed"] and record["removed_area_mm2"] == 4
    assert record["removed_mass_kg"] == pytest.approx(4*5*2.7e-6)
    assert record["geometry_sha256"] != before["geometry_sha256"]
    assert record["rule"]["original_geometry_sha256"] == before["geometry_sha256"]
    assert record["source_field_sha256"] == hashlib.sha256(original.astype("<f8").tobytes()).hexdigest()
    assert record["diagnostics"]["geometry_valid"]
    assert record["diagnostics"]["area_minus_budget_mm2"] == 0
    assert not record["independent_verified"]
    assert record["status"] == "generated_not_analyzed"
    assert record["factorizations"] == record["rhs_solves"] == 0
    assert not repaired.flags.writeable
    json.dumps(record, allow_nan=False)


def test_repair_does_not_invent_pad_connection():
    mesh = make_l_mesh(2)
    with pytest.raises(ValueError, match="both complete pads"):
        retain_pad_component(mesh, bracket_pads(mesh), INPUT_HASH)
    with pytest.raises(ValueError, match="both complete pads"):
        retain_pad_component(mesh, np.zeros(len(mesh.cells)), INPUT_HASH)


@pytest.mark.parametrize("invalid", [[.5]*4, [np.nan]*4, [1, 0], [-1, 0, 0, 1]])
def test_invalid_binary_is_rejected(invalid):
    with pytest.raises(ValueError):
        binary_diagnostics(make_rect_mesh(2, 2, 1), invalid)


@pytest.mark.parametrize("threshold", [-.1, 1.1, np.nan, np.inf])
def test_invalid_threshold_is_rejected(threshold):
    with pytest.raises(ValueError):
        threshold_reconstruction(make_rect_mesh(2, 2, 1), [.1]*4, threshold, INPUT_HASH)


def test_invalid_provenance_and_non_grid_geometry_are_rejected():
    mesh = make_rect_mesh(2, 2, 1)
    with pytest.raises(ValueError, match="SHA-256"):
        source_fingerprint(mesh, np.ones(4), "unverified-display-name")
    distorted = mesh.nodes.copy()
    distorted[4] += [.1, 0]
    irregular = Mesh(distorted, mesh.cells, mesh.cell_indices, 1)
    with pytest.raises(ValueError, match="regular square"):
        binary_diagnostics(irregular, np.ones(4))
    with pytest.raises(ValueError, match="unique"):
        threshold_sweep(mesh, np.ones(4), INPUT_HASH, thresholds=(.5, .5))
