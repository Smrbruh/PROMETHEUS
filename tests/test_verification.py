"""P2 path tests: independent analytical references and smooth convergence."""
from types import SimpleNamespace

import numpy as np
import pytest

from prometheus.models import Material
from prometheus.verification import (
    P2Mesh, assemble_p2, bracket_verify, element_stiffness_p2, error_norms_p2,
    line_force_p2, manufactured_fields, p2_mesh_from_quads, shape_p2, solve_p2,
    triangle_quadrature,
)


def rectangle(width, height, h):
    """Independent test geometry builder; no primary geometry implementation."""
    nx, ny = int(round(width/h)), int(round(height/h))
    nodes = np.array([(i*h, j*h) for j in range(ny+1) for i in range(nx+1)], dtype=float)
    cells = np.array([(j*(nx+1)+i, j*(nx+1)+i+1,
                       (j+1)*(nx+1)+i+1, (j+1)*(nx+1)+i)
                      for j in range(ny) for i in range(nx)])
    return SimpleNamespace(nodes=nodes, cells=cells)


def test_p2_shape_partition_and_nodal_interpolation():
    reference_nodes = [(0, 0), (1, 0), (0, 1), (.5, 0), (.5, .5), (0, .5)]
    for j, (r, s) in enumerate(reference_nodes):
        values, _ = shape_p2(r, s)
        np.testing.assert_allclose(values, np.eye(6)[j], atol=1e-15)
    values, derivatives = shape_p2(.217, .386)
    assert values.sum() == pytest.approx(1.)
    np.testing.assert_allclose(derivatives.sum(axis=0), 0., atol=1e-15)
    for direction in (0, 1):
        point = np.array([.217, .386])
        delta = np.eye(2)[direction]*1e-6
        fd = (shape_p2(*(point+delta))[0]-shape_p2(*(point-delta))[0])/2e-6
        np.testing.assert_allclose(fd, derivatives[:, direction], rtol=1e-9, atol=1e-9)


def test_p2_element_symmetry_three_rigid_modes_and_positive_modes():
    corners = np.array([[0., 0.], [2., .2], [.3, 1.4]])
    points = np.vstack([corners, corners[[[0, 1], [1, 2], [2, 0]]].mean(axis=1)])
    ke = element_stiffness_p2(points, 70000., .3, 5.)
    np.testing.assert_allclose(ke, ke.T, rtol=1e-13, atol=1e-9)
    eigenvalues = np.linalg.eigvalsh(ke)
    scale = eigenvalues[-1]
    assert np.count_nonzero(np.abs(eigenvalues) <= scale*1e-10) == 3
    assert np.count_nonzero(eigenvalues > scale*1e-10) == 9
    modes = np.zeros((12, 3))
    modes[0::2, 0], modes[1::2, 1] = 1, 1
    modes[0::2, 2], modes[1::2, 2] = -points[:, 1], points[:, 0]
    assert np.linalg.norm(ke@modes)/np.linalg.norm(ke) < 1e-14


def test_p2_mesh_assembly_and_midpoint_conformity():
    mesh = p2_mesh_from_quads(rectangle(2, 1, 1))
    assert len(mesh.cells) == 4
    assert len(mesh.nodes) == 15
    assert len(mesh.boundary_edges) == 6
    material = Material()
    moduli = np.array([100., 200.])
    actual = assemble_p2(mesh, moduli, material).toarray()
    expected = np.zeros_like(actual)
    for e, ids in enumerate(mesh.cells):
        dofs = (2*ids[:, None]+[0, 1]).ravel()
        ke = element_stiffness_p2(mesh.nodes[ids], moduli[mesh.parent_cells[e]], material.nu, material.thickness)
        expected[np.ix_(dofs, dofs)] += ke
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-10)


def test_p2_affine_patch_on_irregular_mesh():
    quad = rectangle(3, 3, 1)
    for i, (x, y) in enumerate(quad.nodes):
        if 0 < x < 3 and 0 < y < 3:
            quad.nodes[i] += [.15*np.sin(i), .1*np.cos(i)]
    mesh = p2_mesh_from_quads(quad)
    affine = mesh.nodes@np.array([[.001, -.002], [.003, .0015]])+np.array([.01, -.02])
    boundary = np.unique(mesh.boundary_edges)
    fixed = (2*boundary[:, None]+[0, 1]).ravel()
    result = solve_p2(mesh, 70000., Material(), prescribed=dict(zip(fixed, affine.ravel()[fixed])))
    np.testing.assert_allclose(result.u[:, 0], affine.ravel(), rtol=1e-9, atol=1e-12)
    assert max(result.diagnostics["energy_error_relative"]) < 1e-10


def test_p2_analytical_tensile_strip():
    material = Material()
    mesh = p2_mesh_from_quads(rectangle(100, 10, 5))
    left = np.flatnonzero(mesh.nodes[:, 0] == 0)
    bottom_left = np.flatnonzero(np.all(mesh.nodes == [0, 0], axis=1))[0]
    fixed = {int(2*i): 0. for i in left}
    fixed[int(2*bottom_left+1)] = 0.
    right = mesh.boundary_edges[np.all(mesh.nodes[mesh.boundary_edges, 0] == 100, axis=1)]
    force = line_force_p2(mesh, right, [material.thickness, 0])
    assert force[0::2].sum() == pytest.approx(50.)
    result = solve_p2(mesh, material.E, material, force, fixed)
    analytical = mesh.nodes*np.array([1/material.E, -material.nu/material.E])
    np.testing.assert_allclose(result.u[:, 0], analytical.ravel(), rtol=1e-8, atol=1e-13)
    assert result.compliance[0] == pytest.approx(50*100/material.E, rel=1e-8)
    assert result.diagnostics["reaction_resultant"][0][0] == pytest.approx(-50., rel=1e-10)
    assert max(result.diagnostics["energy_error_relative"]) < 1e-10
    # Every triangle has the same analytical engineering strain and stress.
    from prometheus.verification import _kinematics, _law
    _, grad = shape_p2(1/3, 1/3)
    bmat, _, _ = _kinematics(mesh.nodes[mesh.cells[:, :3]], grad)
    dofs = (2*mesh.cells[:, :, None]+[0, 1]).reshape(-1, 12)
    strain = np.einsum("eij,ej->ei", bmat, result.u[dofs, 0])
    stress = strain@_law(material.E, material.nu).T
    np.testing.assert_allclose(stress, np.tile([1., 0., 0.], (len(mesh.cells), 1)), rtol=1e-8, atol=1e-10)


def test_p2_body_quadrature_area():
    for order in (2, 4, 6):
        assert sum(q[2] for q in triangle_quadrature(order)) == pytest.approx(.5, abs=1e-15)


def test_p2_singular_and_invalid_inputs_fail():
    mesh = p2_mesh_from_quads(rectangle(1, 1, 1))
    material = Material()
    with pytest.raises(ValueError, match="rigid|constraint"):
        solve_p2(mesh, material.E, material)
    with pytest.raises(ValueError, match="rigid|constraint"):
        solve_p2(mesh, material.E, material, prescribed={0: 0., 1: 0.})
    with pytest.raises(ValueError, match="positive"):
        assemble_p2(mesh, 0., material)
    with pytest.raises(ValueError, match="positive"):
        element_stiffness_p2([[0, 0], [0, 1], [1, 0]], material.E, material.nu, 5.)
    with pytest.raises(ValueError, match="incomplete"):
        bracket_verify(rectangle(100, 80, 10), material.E, material)
    # Two squares share a single vertex: an unsupported relative rotation
    # remains although the nodal connectivity graph is connected.
    hinge = SimpleNamespace(
        nodes=np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 1], [2, 2], [1, 2]], dtype=float),
        cells=np.array([[0, 1, 2, 3], [2, 4, 5, 6]]))
    hinge_p2 = p2_mesh_from_quads(hinge)
    with pytest.raises(ValueError, match="singular|mechanism"):
        solve_p2(hinge_p2, material.E, material, prescribed={0: 0., 1: 0., 6: 0., 7: 0.})


def test_p2_independent_bracket_loads_and_superposition():
    # Construct L subset independently of production boundary identification.
    quad = rectangle(120, 80, 10)
    centers = quad.nodes[quad.cells].mean(axis=1)
    cells = quad.cells[(centers[:, 0] < 40) | (centers[:, 1] < 40)]
    used, inverse = np.unique(cells, return_inverse=True)
    quad.nodes, quad.cells = quad.nodes[used], inverse.reshape(-1, 4)
    material = Material()
    result = bracket_verify(quad, material.E, material)
    np.testing.assert_allclose(result.diagnostics["applied_resultant"], [[1000, 0], [0, -1000]], atol=1e-12)
    np.testing.assert_allclose(result.diagnostics["applied_moment"], [-20000, -120000], atol=1e-9)
    assert np.max(np.abs(result.diagnostics["moment_balance"])) < .012
    assert np.max(np.abs(result.diagnostics["resultant_balance"])) < 1e-4
    np.testing.assert_allclose(result.gram, result.gram.T, rtol=1e-10)
    p2 = result.mesh
    top = np.flatnonzero((p2.nodes[:, 1] == 80) & (p2.nodes[:, 0] <= 40))
    prescribed = {int(d): 0. for d in (2*top[:, None]+[0, 1]).ravel()}
    theta = np.deg2rad(13)
    coeff = [np.sin(theta), np.cos(theta)]
    direct = solve_p2(p2, material.E, material, result.forces@coeff, prescribed)
    np.testing.assert_allclose(direct.u[:, 0], result.u@coeff, rtol=1e-8, atol=1e-12)
    assert direct.compliance[0] == pytest.approx(np.array(coeff)@result.gram@coeff, rel=1e-8)


def test_p2_sine_manufactured_solution_convergence_and_quadrature():
    material = Material()
    exact, gradient, force = manufactured_fields(material)
    errors = []
    for h in (5., 2.5, 1.25):
        mesh = p2_mesh_from_quads(rectangle(100, 40, h))
        boundary = np.unique(mesh.boundary_edges)
        prescribed = {int(d): 0. for d in (2*boundary[:, None]+[0, 1]).ravel()}
        result = solve_p2(mesh, material.E, material, prescribed=prescribed, body_force=force)
        e4 = error_norms_p2(mesh, result.u, material, exact, gradient, quadrature_order=4)
        e6 = error_norms_p2(mesh, result.u, material, exact, gradient, quadrature_order=6)
        # Quadrature stability independently of spatial error.
        for key in ("l2", "energy"):
            assert abs(e6[key]-e4[key])/e6[key] < 1e-4
        if h == 5.:
            higher = solve_p2(mesh, material.E, material, prescribed=prescribed,
                              body_force=force, quadrature_order=6)
            assert np.linalg.norm(result.u-higher.u)/np.linalg.norm(higher.u) < 1e-7
        errors.append(e6)
    l2_orders = [np.log2(errors[i]["l2"]/errors[i+1]["l2"]) for i in range(2)]
    energy_orders = [np.log2(errors[i]["energy"]/errors[i+1]["energy"]) for i in range(2)]
    assert min(l2_orders) > 2.85
    assert max(l2_orders) < 3.15
    assert min(energy_orders) > 1.9
    assert max(energy_orders) < 2.1
