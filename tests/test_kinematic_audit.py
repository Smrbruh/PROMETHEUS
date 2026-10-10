"""Affine kinematic identities; no nonlinear equilibrium result is presumed."""
import json
import numpy as np
import pytest

from prometheus.geometry import make_rect_mesh, Mesh
from prometheus.kinematic_audit import audit_q4_kinematics, q4_displacement_gradients


def affine_field(mesh, gradient, translation=(0., 0.)):
    return mesh.nodes @ np.asarray(gradient).T+np.asarray(translation)


def small_mesh():
    return make_rect_mesh(4, 2, 1)


def test_translation_is_exactly_zero_without_zero_strain_ratio():
    mesh = small_mesh()
    u = np.tile([8.125, -12.375], (len(mesh.nodes), 1))
    summary, fields = audit_q4_kinematics(mesh.nodes, mesh.cells, u, singular_points=np.empty((0, 2)))
    assert np.count_nonzero(fields['H']) == 0
    assert np.count_nonzero(fields['green_lagrange']) == 0
    assert summary['global']['quadratic_to_linear_area_l2_ratio'] is None
    assert summary['global']['area_mm2'] == pytest.approx(8)
    assert summary['regular']['area_mm2'] == pytest.approx(8)
    json.dumps(summary, allow_nan=False)


@pytest.mark.parametrize('gradient', [np.array([[1e-4, 0.], [0., -3e-5]]),
                                     np.array([[0., 2e-4], [0., 0.]])])
def test_affine_extension_and_simple_shear(gradient):
    mesh = small_mesh()
    u = affine_field(mesh, gradient)
    summary, fields = audit_q4_kinematics(mesh.nodes, mesh.cells, u)
    np.testing.assert_allclose(fields['H'], np.broadcast_to(gradient, fields['H'].shape), rtol=0, atol=1e-18)
    expected_linear = .5*(gradient+gradient.T)
    expected_quadratic = .5*gradient.T @ gradient
    np.testing.assert_allclose(fields['linear_strain'],
                               np.broadcast_to(expected_linear, fields['linear_strain'].shape), atol=1e-18)
    np.testing.assert_allclose(fields['quadratic_correction'],
                               np.broadcast_to(expected_quadratic, fields['quadratic_correction'].shape), atol=1e-21)
    metric = summary['global']['metrics']['linear_strain_frobenius']
    assert metric['area_l2'] == pytest.approx(np.sqrt(8)*np.linalg.norm(expected_linear), rel=1e-12)
    assert metric['area_rms'] == pytest.approx(np.linalg.norm(expected_linear), rel=1e-12)
    assert metric['area_weighted_p95'] == pytest.approx(np.linalg.norm(expected_linear), rel=1e-12)


@pytest.mark.parametrize('angle', [1e-3, .4, 1.2, -2.4])
def test_exact_finite_rigid_rotation_has_zero_green_lagrange_strain(angle):
    mesh = small_mesh()
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    u = affine_field(mesh, rotation-np.eye(2), translation=(.7, -.3))
    summary, fields = audit_q4_kinematics(mesh.nodes, mesh.cells, u)
    assert np.max(np.abs(fields['green_lagrange'])) < 3e-15
    expected_linear = (np.cos(angle)-1)*np.eye(2)
    np.testing.assert_allclose(fields['linear_strain'],
                               np.broadcast_to(expected_linear, fields['linear_strain'].shape), atol=2e-15)
    assert summary['global']['metrics']['linear_strain_frobenius']['area_rms'] > 0
    assert summary['global']['quadratic_to_linear_area_l2_ratio'] == pytest.approx(1., abs=2e-9)
    np.testing.assert_allclose(fields['polar_rotation_radians'], angle, atol=2e-15)
    np.testing.assert_allclose(fields['infinitesimal_rotation_radians'], np.sin(angle), atol=2e-15)
    np.testing.assert_allclose(fields['determinant_F'], 1., atol=3e-15)


def test_linear_and_quadratic_scaling_are_distinct():
    mesh = small_mesh()
    u = affine_field(mesh, [[.01, .08], [-.03, -.004]])
    first, fields = audit_q4_kinematics(mesh.nodes, mesh.cells, u)
    second, scaled = audit_q4_kinematics(mesh.nodes, mesh.cells, -3*u)
    np.testing.assert_allclose(scaled['H'], -3*fields['H'], atol=3e-16)
    np.testing.assert_allclose(scaled['linear_strain'], -3*fields['linear_strain'], atol=3e-16)
    np.testing.assert_allclose(scaled['quadratic_correction'], 9*fields['quadratic_correction'], atol=3e-17)
    for key in ('area_l2', 'area_rms', 'area_weighted_p95', 'maximum'):
        assert second['global']['metrics']['linear_strain_frobenius'][key] == pytest.approx(
            3*first['global']['metrics']['linear_strain_frobenius'][key])
        assert second['global']['metrics']['quadratic_correction_frobenius'][key] == pytest.approx(
            9*first['global']['metrics']['quadratic_correction_frobenius'][key])
    assert second['global']['quadratic_to_linear_area_l2_ratio'] == pytest.approx(
        3*first['global']['quadratic_to_linear_area_l2_ratio'])


def test_gradient_on_irregular_q4_is_affine_exact():
    nodes = np.array([[0., 0.], [3., .2], [2.8, 1.5], [-.1, 1.]])
    cells = np.array([[0, 1, 2, 3]])
    expected = np.array([[.03, -.02], [.01, .007]])
    u = nodes @ expected.T+np.array([.3, -.2])
    h, points, weights = q4_displacement_gradients(nodes, cells, u, quadrature_order=3)
    np.testing.assert_allclose(h, np.broadcast_to(expected, h.shape), atol=2e-16)
    assert points.shape == (1, 9, 2)
    expected_area = .5*abs(np.sum(nodes[:, 0]*np.roll(nodes[:, 1], -1)-nodes[:, 1]*np.roll(nodes[:, 0], -1)))
    assert weights.sum() == pytest.approx(expected_area, rel=1e-14)


def test_regular_region_is_physical_and_empty_region_is_explicit():
    mesh = make_rect_mesh(4, 2, 1)
    u = affine_field(mesh, [[.01, 0], [0, 0]])
    summary, fields = audit_q4_kinematics(mesh.nodes, mesh.cells, u,
                                         singular_points=np.array([[0., 0.]]), exclusion_radius_mm=2)
    expected = np.linalg.norm(fields['gauss_points_mm'], axis=-1) >= 2
    np.testing.assert_array_equal(fields['regular_region'], expected)
    assert 0 < summary['regular']['area_mm2'] < summary['global']['area_mm2']
    empty, _ = audit_q4_kinematics(mesh.nodes, mesh.cells, u,
                                  singular_points=np.array([[0., 0.]]), exclusion_radius_mm=10)
    assert empty['regular']['status'] == 'not_defined_empty_region'
    assert empty['regular']['metrics'] is None
    json.dumps(empty, allow_nan=False)


def test_reflection_is_flagged_and_has_no_proper_polar_rotation_statistics():
    mesh = small_mesh()
    summary, fields = audit_q4_kinematics(mesh.nodes, mesh.cells,
                                         affine_field(mesh, [[-2., 0.], [0., 0.]]))
    assert summary['global']['nonpositive_determinant_area_fraction'] == 1
    assert summary['global']['polar_rotation_absolute_radians'] is None
    assert np.isnan(fields['polar_rotation_radians']).all()
    json.dumps(summary, allow_nan=False)


def test_area_weighted_percentile_uses_area_not_sample_count():
    # Two disconnected quads are legitimate for a displacement-only diagnostic.
    # A small high-gradient cell has 4% of the area, so P95 stays at the low value.
    nodes = np.array([[0., 0.], [6., 0.], [6., 4.], [0., 4.],
                      [10., 0.], [11., 0.], [11., 1.], [10., 1.]])
    cells = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
    u = np.zeros_like(nodes)
    u[:4, 0] = .01*nodes[:4, 0]
    u[4:, 0] = .1*nodes[4:, 0]
    summary, _ = audit_q4_kinematics(nodes, cells, u)
    metric = summary['global']['metrics']['linear_strain_frobenius']
    assert metric['area_weighted_p95'] == pytest.approx(.01)
    assert metric['area_l2'] == pytest.approx(np.sqrt(24*.01**2+.1**2))


def test_superposed_finite_rotation_preserves_green_lagrange_tensor():
    mesh = small_mesh()
    f = np.array([[1.04, .02], [.01, .98]])
    angle = .6
    r = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    _, before = audit_q4_kinematics(mesh.nodes, mesh.cells, affine_field(mesh, f-np.eye(2)))
    _, after = audit_q4_kinematics(mesh.nodes, mesh.cells, affine_field(mesh, r@f-np.eye(2)))
    np.testing.assert_allclose(after['green_lagrange'], before['green_lagrange'], atol=3e-16)
    assert np.max(np.abs(after['linear_strain']-before['linear_strain'])) > .1


def test_invalid_inputs_are_rejected_and_inputs_are_not_modified():
    mesh = small_mesh()
    u = affine_field(mesh, [[.01, .03], [.02, .01]])
    original = u.copy()
    audit_q4_kinematics(mesh.nodes, mesh.cells, u)
    np.testing.assert_array_equal(u, original)
    with pytest.raises(ValueError, match='one finite displacement'):
        audit_q4_kinematics(mesh.nodes, mesh.cells, np.zeros((mesh.ndof, 2)))
    with pytest.raises(ValueError, match='Quadrature'):
        audit_q4_kinematics(mesh.nodes, mesh.cells, u, quadrature_order=1)
    with pytest.raises(ValueError, match='Exclusion'):
        audit_q4_kinematics(mesh.nodes, mesh.cells, u, exclusion_radius_mm=-1)
    with pytest.raises(ValueError, match='Duplicate'):
        audit_q4_kinematics(mesh.nodes, np.vstack([mesh.cells, mesh.cells[0]]), u)
    with pytest.raises(ValueError, match='Singular points'):
        audit_q4_kinematics(mesh.nodes, mesh.cells, u, singular_points=[np.nan, 0])
