"""Kinematic truncation diagnostics from a saved *linear* Q4 displacement.

For H_ij = du_i/dX_j, the infinitesimal tensor is e=sym(H), the omitted
quadratic term is Q=H.T H/2, and the Green--Lagrange tensor of this same
displacement field is E=e+Q. E is NOT a nonlinear equilibrium solution. The
norm of Q is neither a displacement-error estimate nor a safe-load criterion.

All tensors below are the two-dimensional in-plane block. No out-of-plane
stretch/strain or nonlinear plane-stress constitutive reduction is inferred.
Tensor shear e_xy is used, not engineering shear gamma_xy=2 e_xy.

Primary mechanics references accessed 2026-10-11:
https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.009.html
https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.008.html
https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/resources/lecture-16/

These references distinguish objective finite-deformation strain measures from
the small-displacement/small-rotation strain approximation. The implementation
uses elementary tensor identities, not copied reference code.
"""
import numpy as np

from .fem import quadrature_data, shape_functions


CANONICAL_SINGULAR_POINTS = np.array([[40., 40.], [0., 80.], [40., 80.],
                                    [120., 10.], [120., 30.]])
CANONICAL_SINGULAR_POINTS.setflags(write=False)


def _validate_mesh_and_displacement(nodes, cells, displacement):
    nodes = np.asarray(nodes, dtype=float)
    cells = np.asarray(cells)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or not np.isfinite(nodes).all():
        raise ValueError('Nodes must have finite shape (nnode,2)')
    if (cells.ndim != 2 or cells.shape[1] != 4 or len(cells) == 0
            or not np.issubdtype(cells.dtype, np.integer)):
        raise ValueError('Cells must be nonempty integer Q4 connectivity')
    if np.any(cells < 0) or np.any(cells >= len(nodes)):
        raise ValueError('Cell node index out of range')
    sorted_cells = np.sort(cells, axis=1)
    if np.any(np.diff(sorted_cells, axis=1) == 0):
        raise ValueError('A Q4 cell contains repeated nodes')
    if len(np.unique(sorted_cells, axis=0)) != len(cells):
        raise ValueError('Duplicate cells would double-count reference area')
    displacement = np.asarray(displacement, dtype=float)
    if displacement.shape in ((2*len(nodes),), (2*len(nodes), 1)):
        displacement = displacement.reshape(len(nodes), 2)
    if displacement.shape != nodes.shape or not np.isfinite(displacement).all():
        raise ValueError('Select one finite displacement field: (nnode,2) or (ndof,)')
    return nodes, cells, displacement


def q4_displacement_gradients(nodes, cells, displacement, quadrature_order=2):
    """Return (H, reference points, area weights) with element/point axes.

    Existing Q4 quadrature supplies shape values, reference integration weights
    and Jacobian validity checks. Gradients are independently contracted as
    H=u_nodes.T @ dN/dX, without extracting strain from the solver's B matrix.
    Translationally identical cells share shape-gradient work only.
    """
    nodes, cells, displacement = _validate_mesh_and_displacement(nodes, cells, displacement)
    if (isinstance(quadrature_order, bool) or not isinstance(quadrature_order, (int, np.integer))
            or quadrature_order < 2):
        raise ValueError('Quadrature order must be an integer >=2')
    coordinates = nodes[cells]
    relative = coordinates-coordinates[:, :1]
    shapes, inverse = np.unique(relative.reshape(-1, 8), axis=0, return_inverse=True)
    points_1d, _ = np.polynomial.legendre.leggauss(quadrature_order)
    reference_points = [(xi, eta) for xi in points_1d for eta in points_1d]
    nq = len(reference_points)
    gradients = np.empty((len(shapes), nq, 4, 2))
    weights = np.empty((len(shapes), nq))
    locations = np.empty((len(shapes), nq, 2))
    for shape_index, flattened in enumerate(shapes):
        coords = flattened.reshape(4, 2)
        for point_index, ((xi, eta), (_, _, weight, point)) in enumerate(
                zip(reference_points, quadrature_data(coords, order=quadrature_order))):
            _, reference_gradient = shape_functions(xi, eta)
            jacobian = coords.T @ reference_gradient
            gradients[shape_index, point_index] = reference_gradient @ np.linalg.inv(jacobian)
            weights[shape_index, point_index] = weight
            locations[shape_index, point_index] = point
    # Subtract one element-node displacement: constant translation must cancel
    # exactly instead of leaving a residual from sums of shape gradients.
    local_u = displacement[cells]-displacement[cells[:, :1]]
    h = np.einsum('eai,eqaj->eqij', local_u, gradients[inverse], optimize=True)
    points = locations[inverse]+coordinates[:, :1]
    if not np.isfinite(h).all():
        raise ValueError('Nonfinite displacement gradient')
    return h, points, weights[inverse]


def _weighted_statistics(values, weights):
    """Area-weighted L2 norm, RMS and inverse-CDF percentile of nonnegative scalars."""
    values, weights = np.asarray(values).ravel(), np.asarray(weights).ravel()
    area = float(weights.sum())
    if len(values) == 0 or area == 0:
        return None
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('Statistics require finite values and positive area weights')
    order = np.argsort(values, kind='stable')
    cumulative = np.cumsum(weights[order])
    index = min(int(np.searchsorted(cumulative, .95*area, side='left')), len(values)-1)
    l2 = float(np.sqrt(np.dot(weights, values*values)))
    return {'area_l2': l2, 'area_rms': l2/np.sqrt(area),
            'area_weighted_p95': float(values[order[index]]),
            'maximum': float(values.max())}


def _region_summary(mask, weights, scalar_fields, determinant, polar_angle):
    selected_weights = weights[mask]
    area = float(selected_weights.sum())
    if not np.any(mask):
        return {'status': 'not_defined_empty_region', 'area_mm2': 0., 'gauss_points': 0,
                'metrics': None, 'quadratic_to_linear_area_l2_ratio': None,
                'quadratic_to_linear_p95_ratio': None,
                'ratio_note': 'Empty reference region'}
    metrics = {name: _weighted_statistics(values[mask], selected_weights)
               for name, values in scalar_fields.items()}
    linear = metrics['linear_strain_frobenius']
    correction = metrics['quadratic_correction_frobenius']
    ratio_l2 = correction['area_l2']/linear['area_l2'] if linear['area_l2'] > 0 else None
    ratio_p95 = (correction['area_weighted_p95']/linear['area_weighted_p95']
                 if linear['area_weighted_p95'] > 0 else None)
    polar_defined = mask & (determinant > 0)
    polar_stats = _weighted_statistics(np.abs(polar_angle[polar_defined]), weights[polar_defined])
    orientation_failure_area = float(weights[mask & (determinant <= 0)].sum())
    return {'status': 'diagnostic_only', 'area_mm2': area, 'gauss_points': int(mask.sum()),
            'metrics': metrics, 'quadratic_to_linear_area_l2_ratio': ratio_l2,
            'quadratic_to_linear_p95_ratio': ratio_p95,
            'ratio_note': 'Ratios of region norms/statistics, not a percentile of pointwise ratios; zero denominators return null',
            'polar_rotation_absolute_radians': polar_stats,
            'polar_rotation_defined_area_mm2': float(weights[polar_defined].sum()),
            'determinant_F_minimum': float(determinant[mask].min()),
            'determinant_F_maximum': float(determinant[mask].max()),
            'nonpositive_determinant_area_mm2': orientation_failure_area,
            'nonpositive_determinant_area_fraction': orientation_failure_area/area}


def audit_q4_kinematics(nodes, cells, displacement, quadrature_order=2,
                       exclusion_radius_mm=6., singular_points=None):
    """Return (JSON-ready summary, NPZ-ready fields) for one saved displacement.

    Cells must represent the intended material region; for binary geometries
    use the material-only saved mesh, not an ersatz-filled domain. All area
    weights refer to the original configuration. Global and regular-region
    results are both retained, with no pass/fail load threshold.

    The default regular region is at least 6 mm from the five original v0
    singular points. Pixel reconstructions can introduce additional corners;
    this mask neither identifies nor excludes all of their singularities.
    """
    if not np.isfinite(exclusion_radius_mm) or exclusion_radius_mm < 0:
        raise ValueError('Exclusion radius must be finite and nonnegative')
    singular = CANONICAL_SINGULAR_POINTS if singular_points is None else np.asarray(singular_points, dtype=float)
    if singular.ndim != 2 or singular.shape[1] != 2 or not np.isfinite(singular).all():
        raise ValueError('Singular points must have finite shape (n,2), including (0,2)')
    h, points, weights = q4_displacement_gradients(nodes, cells, displacement, quadrature_order)
    transpose = np.swapaxes(h, -1, -2)
    epsilon = .5*(h+transpose)
    skew = .5*(h-transpose)
    quadratic = .5*np.einsum('eqki,eqkj->eqij', h, h, optimize=True)
    green_lagrange = epsilon+quadratic
    f = h+np.eye(2)
    determinant = np.linalg.det(f)
    if not all(np.isfinite(array).all() for array in (epsilon, skew, quadratic, green_lagrange, determinant)):
        raise ValueError('Nonfinite tensor diagnostic; displacement magnitude may overflow')
    # In 2D this is the proper-polar angle for orientation-preserving F. It is
    # NOT the skew-H angle at finite rotation. Undefined values remain NaN only
    # in NPZ arrays; JSON statistics omit them with explicit area accounting.
    polar_angle = np.full(weights.shape, np.nan)
    preserving = determinant > 0
    polar_angle[preserving] = np.arctan2((f[..., 1, 0]-f[..., 0, 1])[preserving],
                                        (f[..., 0, 0]+f[..., 1, 1])[preserving])
    infinitesimal_rotation = .5*(h[..., 1, 0]-h[..., 0, 1])
    norms = {'displacement_gradient_frobenius': np.linalg.norm(h, axis=(-2, -1)),
             'linear_strain_frobenius': np.linalg.norm(epsilon, axis=(-2, -1)),
             'quadratic_correction_frobenius': np.linalg.norm(quadratic, axis=(-2, -1)),
             'green_lagrange_frobenius': np.linalg.norm(green_lagrange, axis=(-2, -1)),
             'skew_gradient_frobenius': np.linalg.norm(skew, axis=(-2, -1)),
             'infinitesimal_rotation_absolute': np.abs(infinitesimal_rotation)}
    if len(singular):
        regular = np.all(np.linalg.norm(points[:, :, None]-singular[None, None], axis=-1)
                         >= exclusion_radius_mm, axis=-1)
    else:
        regular = np.ones(weights.shape, dtype=bool)
    summary = {
        'schema_version': 1, 'status': 'kinematic_diagnostic_only',
        'input_model': 'One displacement field from the linear model; no new equilibrium solve',
        'tensor_scope': '2D in-plane block only; no thickness strain or nonlinear plane-stress closure',
        'definitions': {'H_ij': 'du_i/dX_j', 'linear_strain': '(H+H.T)/2',
                        'quadratic_correction': 'H.T@H/2',
                        'green_lagrange': 'linear_strain+quadratic_correction',
                        'infinitesimal_rotation': '(H_yx-H_xy)/2; not exact finite rotation',
                        'polar_rotation': 'angle of proper polar rotation of F=I+H when det(F)>0'},
        'statistics': {'weight': 'reference area at Gauss points in mm^2',
                       'area_l2': 'sqrt(sum(w*Frobenius_norm^2)); units mm for dimensionless tensor',
                       'area_rms': 'area_l2/sqrt(reference area); dimensionless',
                       'area_weighted_p95': 'inverse weighted empirical CDF; dimensionless',
                       'rotation_units': 'radians; area_l2 carries radians times mm'},
        'quadrature_order': int(quadrature_order), 'cells': len(cells),
        'exclusion_radius_mm': float(exclusion_radius_mm), 'singular_points_mm': singular.tolist(),
        'regular_region_caveat': 'Original five singularities only; new pixel corners remain and no safety threshold follows',
        'global': _region_summary(np.ones(weights.shape, bool), weights, norms, determinant, polar_angle),
        'regular': _region_summary(regular, weights, norms, determinant, polar_angle),
        'interpretation': 'Quadratic kinematic terms evaluated on a linear field; not nonlinear response, response error, strength, or a safe-load extrapolation',
        'physical_validation': 'not_assessed', 'nonlinear_equilibrium': 'not_solved',
    }
    arrays = {'gauss_points_mm': points, 'reference_area_weights_mm2': weights,
              'H': h, 'linear_strain': epsilon, 'quadratic_correction': quadratic,
              'green_lagrange': green_lagrange, 'skew_gradient': skew,
              'F': f, 'determinant_F': determinant,
              'infinitesimal_rotation_radians': infinitesimal_rotation,
              'polar_rotation_radians': polar_angle, 'regular_region': regular,
              **norms}
    return summary, arrays
