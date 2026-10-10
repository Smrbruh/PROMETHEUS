"""Full-integration Q4 plane-stress elasticity, implemented from weak form.

All arrays use interleaved x/y nodal degrees of freedom. Loads use N,
coordinates mm, moduli MPa, and compliance N mm. Sparse LU is factored once
per solve, with all load columns solved together. No stiffness symmetrization,
regularization, or support insertion is performed.
"""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.sparse.csgraph import connected_components

from .models import NumericalError


@dataclass(frozen=True)
class SolveResult:
    u: np.ndarray
    compliance: np.ndarray
    reactions: np.ndarray
    diagnostics: dict
    wall_seconds: float
    factorizations: int
    rhs_solves: int


def plane_stress(E, nu):
    if not np.isfinite([E, nu]).all() or E <= 0 or not -1 < nu < .5:
        raise ValueError("Elastic modulus must be positive and -1 < nu < 0.5")
    return E/(1-nu*nu)*np.array([[1., nu, 0.], [nu, 1., 0.],
                               [0., 0., (1-nu)/2]])


def shape_functions(xi, eta):
    """Reference Q4 shape functions and derivatives [d/dxi,d/deta]."""
    sx = np.array([-1., 1., 1., -1.])
    sy = np.array([-1., -1., 1., 1.])
    values = .25*(1+sx*xi)*(1+sy*eta)
    deriv = np.column_stack([.25*sx*(1+sy*eta), .25*sy*(1+sx*xi)])
    return values, deriv


def quadrature_data(coords, order=2):
    """Yield (N, B, detJ*weight, physical_point) on a valid Q4 element."""
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (4, 2) or not np.isfinite(coords).all():
        raise ValueError("Q4 coordinates must have finite shape (4,2)")
    points, weights = np.polynomial.legendre.leggauss(order)
    for i, xi in enumerate(points):
        for j, eta in enumerate(points):
            N, dN = shape_functions(xi, eta)
            jac = coords.T @ dN
            det = np.linalg.det(jac)
            scale = max(float(np.linalg.norm(jac))**2, np.finfo(float).tiny)
            if not np.isfinite(det) or det <= 1e-14*scale:
                raise NumericalError("Q4 has nonpositive or degenerate Jacobian", "invalid_geometry")
            gradients = dN @ np.linalg.inv(jac)
            B = np.zeros((3, 8))
            B[0, 0::2], B[1, 1::2] = gradients[:, 0], gradients[:, 1]
            B[2, 0::2], B[2, 1::2] = gradients[:, 1], gradients[:, 0]
            yield N, B, det*weights[i]*weights[j], N @ coords


gauss_data = quadrature_data


def q4_stiffness(coords, E=1., nu=.3, thickness=5.):
    if not np.isfinite(thickness) or thickness <= 0:
        raise ValueError("Thickness must be positive and finite")
    constitutive = plane_stress(E, nu)
    stiffness = np.zeros((8, 8))
    for _, B, weight, _ in quadrature_data(coords):
        stiffness += thickness*weight*(B.T @ constitutive @ B)
    return stiffness


def element_dofs(mesh):
    return (2*mesh.cells[:, :, None]+np.array([0, 1])).reshape(-1, 8)


def unit_stiffness(mesh, nu=.3, thickness=5.):
    """Unit-modulus matrices; reuse exact translationally identical cells."""
    coordinates = mesh.nodes[mesh.cells]
    relative = coordinates-coordinates[:, :1, :]
    shapes, inverse = np.unique(relative.reshape(-1, 8), axis=0, return_inverse=True)
    matrices = np.array([q4_stiffness(shape.reshape(4, 2), 1., nu, thickness)
                         for shape in shapes])
    return matrices[inverse]


def _moduli(mesh, moduli):
    values = np.asarray(moduli, dtype=float)
    if values.ndim == 0:
        values = np.full(len(mesh.cells), values.item())
    if values.shape != (len(mesh.cells),) or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("Supply one strictly positive finite modulus per cell")
    return values


def assemble(mesh, moduli, nu=.3, thickness=5.):
    moduli = _moduli(mesh, moduli)
    edofs = element_dofs(mesh)
    matrices = unit_stiffness(mesh, nu, thickness)*moduli[:, None, None]
    row = np.repeat(edofs, 8, axis=1).ravel()
    col = np.tile(edofs, (1, 8)).ravel()
    K = sparse.coo_matrix((matrices.ravel(), (row, col)),
                          shape=(mesh.ndof, mesh.ndof)).tocsc()
    if not np.isfinite(K.data).all():
        raise NumericalError("Nonfinite assembled stiffness", "invalid_stiffness")
    return K


def _check_support(mesh, fixed):
    """Each node-connected material component must restrain 3 rigid modes."""
    edges = mesh.cells[:, [[0, 1], [1, 2], [2, 3], [3, 0]]].reshape(-1, 2)
    graph = sparse.coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                              shape=(len(mesh.nodes), len(mesh.nodes))).tocsr()
    count, labels = connected_components(graph, directed=False)
    for group in range(count):
        nodeids = np.flatnonzero(labels == group)
        constrained = fixed[labels[fixed//2] == group]
        if len(constrained) < 3:
            raise NumericalError("Material component has unconstrained rigid motion", "singular_system")
        center = mesh.nodes[nodeids].mean(axis=0)
        xy = mesh.nodes[constrained//2]-center
        length = max(np.ptp(mesh.nodes[nodeids], axis=0).max(), 1.)
        modes = np.zeros((len(constrained), 3))
        modes[:, 0] = constrained % 2 == 0
        modes[:, 1] = constrained % 2 == 1
        modes[:, 2] = np.where(constrained % 2 == 0, -xy[:, 1], xy[:, 0])/length
        if np.linalg.matrix_rank(modes, tol=1e-12) != 3:
            raise NumericalError("Support leaves a rigid translation or rotation", "singular_system")


def solve(mesh, moduli, forces, fixed_dofs, nu=.3, thickness=5., prescribed=None):
    start = perf_counter()
    forces = np.asarray(forces, dtype=float)
    if forces.ndim == 1:
        forces = forces[:, None]
    if (forces.ndim != 2 or forces.shape[0] != mesh.ndof or forces.shape[1] == 0
            or not np.isfinite(forces).all()):
        raise ValueError("Forces must be finite with shape (ndof,nrhs) or (ndof,)")
    raw_fixed = np.asarray(fixed_dofs)
    if raw_fixed.ndim != 1 or not np.issubdtype(raw_fixed.dtype, np.integer):
        raise ValueError("Constrained DOFs must be a one-dimensional integer array")
    fixed = np.unique(raw_fixed.astype(np.int64))
    if len(fixed) == 0:
        raise NumericalError("A free body has no supported solution", "singular_system")
    if np.any(fixed < 0) or np.any(fixed >= mesh.ndof):
        raise ValueError("Constrained degree of freedom is out of range")
    _check_support(mesh, fixed)
    free = np.setdiff1d(np.arange(mesh.ndof), fixed)
    u = np.zeros_like(forces)
    if prescribed is not None:
        values = np.asarray(prescribed, dtype=float)
        if values.ndim == 1:
            values = values[:, None]
        if values.ndim != 2:
            raise ValueError("Prescribed displacements must be a vector or matrix")
        if values.shape[0] == mesh.ndof:
            values = values[fixed]
        if values.shape[0] != len(fixed) or values.shape[1] not in (1, forces.shape[1]):
            raise ValueError("Prescribed values must match fixed or all DOFs and RHS count")
        if not np.isfinite(values).all():
            raise ValueError("Prescribed displacement must be finite")
        u[fixed] = values
    K = assemble(mesh, moduli, nu, thickness)
    norm_k = np.linalg.norm(K.data)
    symmetry_error = np.linalg.norm((K-K.T).data)/max(norm_k, 1.)
    rhs = forces[free]-K[free][:, fixed] @ u[fixed]
    factorizations = 0
    pivot_ratio = None
    if len(free):
        try:
            factor = splu(K[free][:, free].tocsc(), permc_spec="COLAMD")
            pivots = np.abs(factor.U.diagonal())
            pivot_ratio = float(pivots.min()/pivots.max())
            if not np.isfinite(pivot_ratio) or pivot_ratio < 1e-14:
                raise NumericalError("Stiffness is singular or numerically near-singular", "singular_system")
            u[free] = factor.solve(rhs)
            factorizations = 1
        except RuntimeError as exc:
            raise NumericalError("Sparse factorization failed: "+str(exc), "singular_system") from exc
    if not np.isfinite(u).all():
        raise NumericalError("Nonfinite solution", "nonfinite_solution")
    ku = K @ u
    full_residual = ku-forces
    residual = np.linalg.norm(full_residual[free], axis=0)/np.maximum(np.linalg.norm(rhs, axis=0), 1.)
    reactions = np.zeros_like(forces)
    reactions[fixed] = full_residual[fixed]
    compliance = np.einsum("ij,ij->j", forces, u)
    strain_energy_twice = np.einsum("ij,ij->j", u, ku)
    # Prescribed displacement work is essential for nonzero Dirichlet data.
    total_work = compliance+np.einsum("ij,ij->j", reactions, u)
    energy_error = np.abs(total_work-strain_energy_twice)/np.maximum(np.abs(strain_energy_twice), 1e-12)
    element_u = u[element_dofs(mesh)]
    element_energy = np.einsum("eir,eij,ejr,e->r", element_u,
                              unit_stiffness(mesh, nu, thickness), element_u,
                              _moduli(mesh, moduli), optimize=True)
    element_energy_error = np.abs(element_energy-strain_energy_twice)/np.maximum(np.abs(strain_energy_twice), 1e-12)
    nrhs = forces.shape[1]
    applied = forces.reshape(-1, 2, nrhs).sum(axis=0)
    reaction_resultant = reactions.reshape(-1, 2, nrhs).sum(axis=0)
    scale = np.maximum(np.maximum(np.linalg.norm(applied, axis=0),
                                  np.linalg.norm(reaction_resultant, axis=0)), 1.)
    balance = np.linalg.norm(applied+reaction_resultant, axis=0)/scale
    total = (forces+reactions).reshape(-1, 2, nrhs)
    moment = np.sum(mesh.nodes[:, 0, None]*total[:, 1]
                    - mesh.nodes[:, 1, None]*total[:, 0], axis=0)
    moment_scale = np.maximum(scale*max(np.ptp(mesh.nodes, axis=0).max(), 1.), 1.)
    moment_error = np.abs(moment)/moment_scale
    homogeneous = not np.any(u[fixed])
    positive = np.all(compliance >= 0) if homogeneous else True
    diagnostics = {
        "residual_max": float(residual.max()),
        "symmetry_error": float(symmetry_error),
        "energy_error_max": float(energy_error.max()),
        "element_energy_error_max": float(element_energy_error.max()),
        "force_balance_max": float(balance.max()),
        "moment_balance_max": float(moment_error.max()),
        "minimum_relative_lu_pivot": pivot_ratio,
        "pivot_ratio_is_condition_estimate": False,
        "energy_includes_prescribed_displacement_work": not homogeneous,
        "solver": "scipy.sparse.linalg.splu", "ordering": "COLAMD",
        "passed": bool(np.all(residual <= 1e-8) and symmetry_error <= 1e-12
                        and np.all(energy_error <= 1e-8) and np.all(balance <= 1e-7)
                        and np.all(element_energy_error <= 1e-8)
                        and np.all(moment_error <= 1e-7) and positive),
    }
    if not diagnostics["passed"]:
        raise NumericalError("Solution failed equilibrium/energy diagnostics: "+str(diagnostics),
                             "verification_failed")
    return SolveResult(u, compliance, reactions, diagnostics, perf_counter()-start,
                       factorizations, nrhs if len(free) else 0)


def stress_metrics(mesh, moduli, u, nu=.3, thickness=5., densities=None):
    """Diagnostic integration-point von Mises stresses (no strength criterion).

    The regular region is at least 6 mm from the five specified singular points
    and in cells with density >= 0.9. Its percentile is weighted by physical area.
    Uses a single displacement vector; select a load case before calling.
    """
    moduli = _moduli(mesh, moduli)
    u = np.asarray(u, dtype=float).reshape(-1)
    if len(u) != mesh.ndof or not np.isfinite(u).all():
        raise ValueError("Stress requires one finite displacement vector")
    if densities is None:
        densities = np.ones(len(mesh.cells))
    densities = np.asarray(densities, dtype=float)
    if densities.shape != (len(mesh.cells),) or not np.isfinite(densities).all():
        raise ValueError("Densities must have shape (ne,)")
    edofs = element_dofs(mesh)
    singular = np.array([[40., 40.], [0., 80.], [40., 80.], [120., 10.], [120., 30.]])
    vm = np.empty((len(mesh.cells), 4))
    positions = np.empty((len(mesh.cells), 4, 2))
    areas = np.empty_like(vm)
    D = plane_stress(1., nu)
    coordinates = mesh.nodes[mesh.cells]
    relative = coordinates-coordinates[:, :1, :]
    shapes, inverse = np.unique(relative.reshape(-1, 8), axis=0, return_inverse=True)
    local_u = u[edofs]
    for shape_id, shape in enumerate(shapes):
        selected = np.flatnonzero(inverse == shape_id)
        for g, (N, B, weight, _) in enumerate(quadrature_data(shape.reshape(4, 2))):
            stress = (local_u[selected] @ B.T @ D.T)*moduli[selected, None]
            sx, sy, tau = stress.T
            vm[selected, g] = np.sqrt(np.maximum(sx*sx-sx*sy+sy*sy+3*tau*tau, 0.))
            positions[selected, g] = np.einsum("i,eij->ej", N, coordinates[selected])
            areas[selected, g] = weight
    distance = np.linalg.norm(positions[:, :, None, :]-singular[None, None, :, :], axis=3).min(axis=2)
    regular = (distance >= 6.) & (densities[:, None] >= .9)
    material = densities > 0
    peak = float(vm[material].max()) if material.any() else None
    if regular.any():
        values, weights = vm[regular], areas[regular]
        order = np.argsort(values, kind="stable")
        cumulative = np.cumsum(weights[order])
        percentile = float(values[order[np.searchsorted(cumulative, .95*cumulative[-1])]])
    else:
        percentile = None
    return {"unrestricted_peak_mpa": peak,
            "regular_p95_mpa": percentile,
            "regular_region_area_mm2": float(areas[regular].sum()),
            "regular_region_status": "defined" if regular.any() else "not_defined",
            "von_mises_gauss_mpa": vm, "gauss_coordinates_mm": positions,
            "diagnostic_only": True}
