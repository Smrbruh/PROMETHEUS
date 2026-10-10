"""Analytic, geometric, and mechanical invariants for the primary Q4 solver."""
import numpy as np
import pytest

from prometheus.geometry import (Mesh, make_rect_mesh, make_l_mesh, subset_mesh,
    bracket_pads, bracket_load_basis, bracket_fixed_dofs)
from prometheus.fem import (plane_stress, q4_stiffness, assemble, solve,
    quadrature_data, element_dofs, stress_metrics)
from prometheus.models import NumericalError


def test_canonical_geometry_and_integrated_patch_load():
    mesh = make_l_mesh()
    assert len(mesh.cells) == 1600
    assert mesh.areas.sum() == 6400
    assert np.all(mesh.areas == 4)
    assert bracket_pads(mesh).sum() == 150
    assert len(np.unique(mesh.cells)) == len(mesh.nodes)
    fixed = bracket_fixed_dofs(mesh)
    assert len(fixed) == 42
    assert np.all(mesh.nodes[fixed//2, 1] == 80)
    B = bracket_load_basis(mesh)
    np.testing.assert_allclose(B.reshape(-1, 2, 2).sum(axis=0), [[1, 0], [0, -1]], atol=1e-15)
    moments = (mesh.nodes[:, 0, None]*B[1::2]-mesh.nodes[:, 1, None]*B[0::2]).sum(axis=0)
    np.testing.assert_allclose(moments, [-20, -120], atol=1e-13)
    endpoint = np.flatnonzero(np.all(mesh.nodes == [120, 10], axis=1))[0]
    interior = np.flatnonzero(np.all(mesh.nodes == [120, 12], axis=1))[0]
    assert B[2*endpoint, 0] == pytest.approx(.05)
    assert B[2*interior, 0] == pytest.approx(.1)


def test_load_and_support_patch_must_be_complete_actual_boundary():
    mesh = make_l_mesh()
    keep = ~((mesh.centroids[:, 0] == 119) & (mesh.centroids[:, 1] == 11))
    with pytest.raises(ValueError, match="patch"):
        bracket_load_basis(subset_mesh(mesh, keep))
    keep = ~((mesh.centroids[:, 0] == 1) & (mesh.centroids[:, 1] == 79))
    with pytest.raises(ValueError, match="patch"):
        bracket_fixed_dofs(subset_mesh(mesh, keep))


def test_element_exact_integrated_square_and_rigid_modes():
    coords = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    actual = q4_stiffness(coords, E=1, nu=0, thickness=1)
    # Exact separable integrals: int_0^1 x^2=1/3, int x(1-x)=1/6.
    xx = np.array([[2, -1, -1, 0], [-1, 2, 0, -1],
                   [-1, 0, 2, -1], [0, -1, -1, 2]])/4
    yy = np.array([[2, 0, -1, -1], [0, 2, -1, -1],
                   [-1, -1, 2, 0], [-1, -1, 0, 2]])/4
    xy = np.array([[1, -1, -1, 1], [1, -1, -1, 1],
                   [-1, 1, 1, -1], [-1, 1, 1, -1]])/8
    exact = np.empty((8, 8))
    exact[0::2, 0::2], exact[1::2, 1::2] = xx, yy
    exact[0::2, 1::2], exact[1::2, 0::2] = xy, xy.T
    np.testing.assert_allclose(actual, exact, atol=2e-15)
    np.testing.assert_allclose(actual, actual.T, atol=1e-15)
    eigenvalues = np.linalg.eigvalsh(actual)
    assert np.max(np.abs(eigenvalues[:3])) < 1e-10*eigenvalues[-1]
    assert np.all(eigenvalues[3:] > 1e-10*eigenvalues[-1])
    rigid = np.column_stack([np.tile([1, 0], 4), np.tile([0, 1], 4),
                             np.column_stack([-coords[:, 1], coords[:, 0]]).ravel()])
    np.testing.assert_allclose(actual @ rigid, 0, atol=2e-15)


def test_assembly_matches_explicit_element_scatter():
    mesh = make_rect_mesh(2, 2, 1)
    moduli = np.arange(1., 5.)
    actual = assemble(mesh, moduli).toarray()
    expected = np.zeros_like(actual)
    for cell, ids, E in zip(mesh.cells, element_dofs(mesh), moduli):
        expected[np.ix_(ids, ids)] += q4_stiffness(mesh.nodes[cell], E)
    np.testing.assert_allclose(actual, expected, atol=1e-14)
    np.testing.assert_allclose(actual, actual.T, atol=1e-14)


def test_irregular_constant_strain_patch_nonzero_boundary_work():
    regular = make_rect_mesh(2, 2, 1)
    nodes = regular.nodes.copy()
    nodes[4] += [.21, -.16]
    mesh = Mesh(nodes, regular.cells, regular.cell_indices, 1.)
    A = np.array([[.012, -.003], [.004, -.007]])
    exact = (nodes @ A.T + np.array([.004, -.006])).ravel()
    boundary = np.flatnonzero((nodes[:, 0] == 0) | (nodes[:, 0] == 2)
                               | (nodes[:, 1] == 0) | (nodes[:, 1] == 2))
    fixed = (2*boundary[:, None]+[0, 1]).ravel()
    result = solve(mesh, 70000., np.zeros(mesh.ndof), fixed, prescribed=exact)
    np.testing.assert_allclose(result.u[:, 0], exact, rtol=1e-9, atol=1e-12)
    expected_strain = [A[0, 0], A[1, 1], A[0, 1]+A[1, 0]]
    for cell, ids in zip(mesh.cells, element_dofs(mesh)):
        for _, B, _, _ in quadrature_data(mesh.nodes[cell]):
            np.testing.assert_allclose(B @ result.u[ids, 0], expected_strain, atol=1e-12)
    assert result.diagnostics["passed"]
    assert result.compliance[0] == 0
    assert result.diagnostics["energy_includes_prescribed_displacement_work"]


def test_analytical_uniaxial_tensile_strip():
    mesh = make_rect_mesh(100, 10, 2)
    right = np.flatnonzero(mesh.nodes[:, 0] == 100)
    force = np.zeros(mesh.ndof)
    force[2*right] = 10.
    force[2*right[[0, -1]]] = 5.
    fixed = np.r_[2*np.flatnonzero(mesh.nodes[:, 0] == 0), 1]
    result = solve(mesh, 70000., force, fixed)
    exact = np.column_stack([mesh.nodes[:, 0]/70000., -.3*mesh.nodes[:, 1]/70000.]).ravel()
    np.testing.assert_allclose(result.u[:, 0], exact, rtol=1e-8, atol=1e-13)
    assert force.sum() == 50
    assert result.reactions[0::2].sum() == pytest.approx(-50, rel=1e-8)
    assert result.compliance[0] == pytest.approx(50*100/70000, rel=1e-8)
    stresses = stress_metrics(mesh, 70000., result.u[:, 0])
    np.testing.assert_allclose(stresses["von_mises_gauss_mpa"], 1., rtol=1e-8)


def test_superposition_scaling_and_energy():
    mesh = make_l_mesh(5)
    B = bracket_load_basis(mesh)*1000.
    fixed = bracket_fixed_dofs(mesh)
    result = solve(mesh, 70000., B, fixed)
    angle = np.deg2rad(13.)
    coefficients = [np.sin(angle), np.cos(angle)]
    direct = solve(mesh, 70000., B @ coefficients, fixed)
    np.testing.assert_allclose(direct.u[:, 0], result.u @ coefficients, rtol=1e-8, atol=1e-12)
    Q = B.T @ result.u
    np.testing.assert_allclose(Q, Q.T, rtol=1e-10)
    assert direct.compliance[0] == pytest.approx(np.array(coefficients) @ Q @ coefficients, rel=1e-8)
    doubled_load = solve(mesh, 70000., 2*B[:, 1], fixed)
    doubled_E = solve(mesh, 140000., B[:, 1], fixed)
    np.testing.assert_allclose(doubled_load.u[:, 0], 2*result.u[:, 1], rtol=1e-9, atol=1e-12)
    assert doubled_load.compliance[0] == pytest.approx(4*result.compliance[1], rel=1e-9)
    np.testing.assert_allclose(doubled_E.u[:, 0], .5*result.u[:, 1], rtol=1e-9, atol=1e-12)
    assert doubled_E.compliance[0] == pytest.approx(.5*result.compliance[1], rel=1e-9)
    assert result.factorizations == 1 and result.rhs_solves == 2
    assert result.diagnostics["passed"]
    assert np.max(np.abs(result.u[fixed])) == 0


def test_singular_and_invalid_input_fail_explicitly():
    mesh = make_rect_mesh(2, 2, 1)
    with pytest.raises(NumericalError, match="free body"):
        solve(mesh, 1., np.zeros(mesh.ndof), np.array([], dtype=int))
    with pytest.raises(NumericalError, match="rigid"):
        solve(mesh, 1., np.zeros(mesh.ndof), np.array([0, 1]))
    with pytest.raises(ValueError, match="modulus"):
        assemble(mesh, 0.)
    with pytest.raises(ValueError, match="finite"):
        solve(mesh, 1., np.full(mesh.ndof, np.nan), np.array([0, 1, 2]))
    with pytest.raises(NumericalError, match="Jacobian"):
        q4_stiffness(mesh.nodes[mesh.cells[0][::-1]])
    with pytest.raises(ValueError):
        plane_stress(1., .5)
    with pytest.raises(ValueError):
        make_rect_mesh(2, 2, .3)
    # A retained disconnected island may not acquire artificial support.
    full = make_rect_mesh(4, 1, 1)
    disconnected = subset_mesh(full, np.array([True, False, False, True]))
    fixed = (2*np.flatnonzero(disconnected.nodes[:, 0] == 0)[:, None]+[0, 1]).ravel()
    with pytest.raises(NumericalError, match="rigid"):
        solve(disconnected, 1., np.zeros(disconnected.ndof), fixed)


def test_empty_regular_stress_region_is_not_defined():
    mesh = make_rect_mesh(2, 2, 1)
    result = stress_metrics(mesh, 1., np.zeros(mesh.ndof), densities=np.full(4, .5))
    assert result["regular_p95_mpa"] is None
    assert result["regular_region_area_mm2"] == 0
    assert result["regular_region_status"] == "not_defined"


def manufactured_errors(h, order=4):
    """Specified smooth fixture, independently analytical body force/derivatives.

    Vectorized over uniform cells to keep this convergence test inexpensive.
    Returns absolute displacement L2 and strain-energy errors, plus solver data.
    """
    E, nu, thickness = 70000., .3, 5.
    a, b, kx, ky = .01, .02, np.pi/100, np.pi/40
    lam, mu = E*nu/(1-nu**2), E/(2*(1+nu))
    mesh = make_rect_mesh(100, 40, h)
    coords = mesh.nodes[mesh.cells]
    edofs = element_dofs(mesh)
    force = np.zeros(mesh.ndof)
    local_force = np.zeros((len(mesh.cells), 8))
    quadrature = list(quadrature_data(coords[0], order))
    for N, _, weight, _ in quadrature:
        xy = np.einsum("i,eij->ej", N, coords)
        x, y = xy.T
        phi = np.sin(kx*x)*np.sin(ky*y)
        chi = np.cos(kx*x)*np.cos(ky*y)
        fx = a*((lam+2*mu)*kx*kx+mu*ky*ky)*phi-(lam+mu)*b*kx*ky*chi
        fy = b*(mu*kx*kx+(lam+2*mu)*ky*ky)*phi-(lam+mu)*a*kx*ky*chi
        local_force[:, 0::2] += thickness*weight*fx[:, None]*N
        local_force[:, 1::2] += thickness*weight*fy[:, None]*N
    np.add.at(force, edofs.ravel(), local_force.ravel())
    x, y = mesh.nodes.T
    boundary = np.flatnonzero((x == 0) | (x == 100) | (y == 0) | (y == 40))
    fixed = (2*boundary[:, None]+[0, 1]).ravel()
    solution = solve(mesh, E, force, fixed)
    local_u = solution.u[edofs, 0]
    l2_squared, energy_squared = 0., 0.
    D = plane_stress(E, nu)
    for N, B, weight, _ in quadrature:
        xy = np.einsum("i,eij->ej", N, coords)
        x, y = xy.T
        phi = np.sin(kx*x)*np.sin(ky*y)
        exact_u = phi[:, None]*[a, b]
        numerical_u = np.einsum("i,eij->ej", N, local_u.reshape(-1, 4, 2))
        dx = kx*np.cos(kx*x)*np.sin(ky*y)
        dy = ky*np.sin(kx*x)*np.cos(ky*y)
        exact_strain = np.column_stack([a*dx, b*dy, a*dy+b*dx])
        numerical_strain = local_u @ B.T
        strain_error = numerical_strain-exact_strain
        l2_squared += weight*np.sum((numerical_u-exact_u)**2)
        energy_squared += thickness*weight*np.einsum("ei,ij,ej->", strain_error, D, strain_error)
    return {"h_mm": h, "l2_error": float(np.sqrt(l2_squared)),
            "energy_error": float(np.sqrt(energy_squared)),
            "quadrature_order": order, "solver": solution.diagnostics}


def test_smooth_manufactured_solution_convergence_and_quadrature_stability():
    measured = [manufactured_errors(h, 4) for h in (5., 2.5, 1.25)]
    errors = np.array([[row["l2_error"], row["energy_error"]] for row in measured])
    orders = np.log2(errors[:-1]/errors[1:])
    # Q4 smooth-solution asymptotic rates are two in L2 and one in energy.
    np.testing.assert_allclose(orders[:, 0], 2., atol=.08)
    np.testing.assert_allclose(orders[:, 1], 1., atol=.04)
    higher = manufactured_errors(1.25, 6)
    np.testing.assert_allclose([higher["l2_error"], higher["energy_error"]],
                                errors[-1], rtol=1e-6)
