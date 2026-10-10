"""Independent straight-sided quadratic-triangle elasticity verifier.

This module shares no quadrilateral shape functions, constitutive matrices,
assembly routines, boundary tags, or load vectors with the primary solver.
It does share SciPy's sparse direct linear algebra. A quadrilateral is split
along corner 0--2 and its frozen modulus is copied to both P2 triangles.
All calculations use mm, N, and plane stress; forces are line loads in N/mm.
"""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu


@dataclass(frozen=True)
class P2Mesh:
    nodes: np.ndarray
    cells: np.ndarray  # corners 0,1,2, followed by midpoints 01,12,20
    parent_cells: np.ndarray
    boundary_edges: np.ndarray  # endpoint, midpoint, endpoint

    @property
    def ndof(self):
        return 2 * len(self.nodes)

    @property
    def ntriangles(self):
        return len(self.cells)


@dataclass
class P2Result:
    mesh: P2Mesh
    u: np.ndarray  # (ndof,nload), even for one right-hand side
    forces: np.ndarray
    compliance: np.ndarray
    gram: np.ndarray
    diagnostics: dict
    wall_seconds: float
    max_pad_displacement: float = float("nan")


def p2_mesh_from_quads(mesh):
    """Promote conforming CCW quads to continuous P2 triangles independently."""
    xy = np.asarray(mesh.nodes, dtype=float)
    quads = np.asarray(mesh.cells)
    if (xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all()
            or quads.ndim != 2 or quads.shape[1] != 4 or len(quads) == 0
            or not np.issubdtype(quads.dtype, np.integer)):
        raise ValueError("Expected finite points and integer four-node cells")
    if np.any(quads < 0) or np.any(quads >= len(xy)):
        raise ValueError("Invalid quadrilateral node index")
    tri = quads[:, [[0, 1, 2], [0, 2, 3]]].reshape(-1, 3)
    corners = xy[tri]
    sides = corners[:, 1:] - corners[:, :1]
    det = sides[:, 0, 0]*sides[:, 1, 1] - sides[:, 0, 1]*sides[:, 1, 0]
    if np.any(det <= 0):
        raise ValueError("P2 triangle Jacobians must be strictly positive")
    all_edges = tri[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2)
    edges, inv, counts = np.unique(np.sort(all_edges, axis=1), axis=0,
                                   return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("Nonmanifold triangle edge")
    mid_ids = np.arange(len(edges)) + len(xy)
    cells = np.column_stack([tri, mid_ids[inv].reshape(-1, 3)])
    boundary = np.column_stack([edges[:, 0], mid_ids, edges[:, 1]])[counts == 1]
    return P2Mesh(np.vstack([xy, xy[edges].mean(axis=1)]), cells,
                  np.repeat(np.arange(len(quads)), 2), boundary)


def shape_p2(r, s):
    """Six Lagrange values and their reference derivatives (6,2)."""
    bary = np.array([1.0-r-s, r, s])
    db = np.array([[-1., -1.], [1., 0.], [0., 1.]])
    vals = np.empty(6)
    deriv = np.empty((6, 2))
    vals[:3] = bary * (2*bary-1)
    deriv[:3] = (4*bary-1)[:, None]*db
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0)), start=3):
        vals[k] = 4*bary[i]*bary[j]
        deriv[k] = 4*(bary[i]*db[j]+bary[j]*db[i])
    return vals, deriv


def _law(modulus, nu):
    """Independent Lamé construction of the engineering-strain law."""
    if not np.isfinite(nu) or not -1 < nu < .5:
        raise ValueError("Poisson ratio must be finite and lie in (-1,.5)")
    lam = modulus*nu/(1-nu*nu)
    mu = modulus/(2*(1+nu))
    return np.array([[lam+2*mu, lam, 0.], [lam, lam+2*mu, 0.], [0., 0., mu]])


def _kinematics(corners, deriv):
    # J columns are physical vectors for the reference r and s directions.
    jac = np.stack([corners[:, 1]-corners[:, 0],
                    corners[:, 2]-corners[:, 0]], axis=2)
    det = np.linalg.det(jac)
    if not np.isfinite(det).all() or np.any(det <= 0):
        raise ValueError("P2 triangle Jacobians must be positive")
    gradient = np.einsum("ai,eij->eaj", deriv, np.linalg.inv(jac))
    bmat = np.zeros((len(corners), 3, 12))
    bmat[:, 0, 0::2] = gradient[:, :, 0]
    bmat[:, 1, 1::2] = gradient[:, :, 1]
    bmat[:, 2, 0::2] = gradient[:, :, 1]
    bmat[:, 2, 1::2] = gradient[:, :, 0]
    return bmat, gradient, det


_STIFFNESS_RULE = ((1/6, 1/6, 1/6), (2/3, 1/6, 1/6), (1/6, 2/3, 1/6))


def element_stiffness_p2(coordinates, modulus, nu, thickness):
    """12x12 stiffness, exactly integrated by a degree-two triangle rule."""
    coordinates = np.asarray(coordinates, dtype=float)
    if (coordinates.shape not in ((3, 2), (6, 2)) or not np.isfinite(coordinates).all()
            or not np.isfinite([modulus, thickness]).all() or min(modulus, thickness) <= 0):
        raise ValueError("Invalid triangle, modulus, or thickness")
    if coordinates.shape == (6, 2):
        mids = coordinates[[[0, 1], [1, 2], [2, 0]]].mean(axis=1)
        if not np.allclose(coordinates[3:], mids, rtol=0, atol=1e-10):
            raise ValueError("Only straight-sided P2 triangles are supported")
    law = _law(modulus, nu)
    result = np.zeros((12, 12))
    for r, s, weight in _STIFFNESS_RULE:
        _, deriv = shape_p2(r, s)
        bmat, _, det = _kinematics(coordinates[None, :3], deriv)
        result += thickness*weight*det[0]*(bmat[0].T @ law @ bmat[0])
    return result


def _triangle_moduli(mesh, moduli):
    values = np.asarray(moduli, dtype=float)
    if values.ndim == 0:
        values = np.full(len(mesh.cells), float(values))
    elif values.shape == (int(mesh.parent_cells.max())+1,):
        values = values[mesh.parent_cells]
    elif values.shape != (len(mesh.cells),):
        raise ValueError("Modulus must be scalar, per parent quad, or per triangle")
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("Moduli must be positive and finite")
    return values


def assemble_p2(mesh, moduli, material):
    """Assemble independent P2 stiffness; parent coefficients remain frozen."""
    moduli = _triangle_moduli(mesh, moduli)
    thickness = float(material.thickness)
    if not np.isfinite(thickness) or thickness <= 0:
        raise ValueError("Thickness must be positive and finite")
    law_unit = _law(1., material.nu)
    # COO assembly itself scales as 144 entries per triangle. This conservative
    # cap avoids initiating an allocation likely to violate the 8 GiB PC budget;
    # it does not predict sparse-factor fill or replace runner resource checks.
    estimate = len(mesh.cells)*144*40 + mesh.ndof*64
    if estimate > 4*1024**3:
        raise MemoryError("P2 assembly estimate exceeds 4 GiB; factorization needs headroom")
    dofs = (2*mesh.cells[:, :, None]+np.array([0, 1])).reshape(-1, 12)
    data = np.empty((len(mesh.cells), 12, 12))
    for start in range(0, len(mesh.cells), 2048):
        end = min(start+2048, len(mesh.cells))
        corners = mesh.nodes[mesh.cells[start:end, :3]]
        block = np.zeros((end-start, 12, 12))
        for r, s, weight in _STIFFNESS_RULE:
            _, deriv = shape_p2(r, s)
            bmat, _, det = _kinematics(corners, deriv)
            block += np.einsum("eai,ab,ebj,e->eij", bmat, law_unit, bmat,
                               moduli[start:end]*det*thickness*weight, optimize=True)
        data[start:end] = block
    rows = np.repeat(dofs, 12, axis=1).ravel()
    cols = np.tile(dofs, (1, 12)).ravel()
    return sparse.coo_matrix((data.ravel(), (rows, cols)),
                             shape=(mesh.ndof, mesh.ndof)).tocsr()


def triangle_quadrature(order=4):
    """Tensor Gauss rule via Duffy map; weights integrate reference area 1/2."""
    if not isinstance(order, (int, np.integer)) or order < 2:
        raise ValueError("Triangle quadrature order must be an integer >=2")
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x+1)/2, w/2
    return [(r, (1-r)*v, wr*wv*(1-r)) for r, wr in zip(x, w) for v, wv in zip(x, w)]


def body_force_p2(mesh, material, body_force, quadrature_order=4):
    """Integrate a physical body force in N/mm^3, including plate thickness."""
    force = np.zeros(mesh.ndof)
    corners = mesh.nodes[mesh.cells[:, :3]]
    for r, s, weight in triangle_quadrature(quadrature_order):
        vals, deriv = shape_p2(r, s)
        _, _, det = _kinematics(corners, deriv)
        xy = (1-r-s)*corners[:, 0]+r*corners[:, 1]+s*corners[:, 2]
        physical = np.asarray(body_force(xy), dtype=float)
        if physical.shape != (len(corners), 2) or not np.isfinite(physical).all():
            raise ValueError("Body force must return a finite (n,2) array")
        terms = (physical[:, None, :]*vals[None, :, None]
                 *(det*weight*material.thickness)[:, None, None])
        np.add.at(force, (2*mesh.cells[:, :, None]+[0, 1]).ravel(), terms.ravel())
    return force


def line_force_p2(mesh, edges, traction):
    """Exact consistent force for uniform N/mm traction on quadratic edges."""
    edges = np.asarray(edges)
    traction = np.asarray(traction, dtype=float)
    if edges.ndim != 2 or edges.shape[1] != 3 or traction.shape != (2,) or not np.isfinite(traction).all():
        raise ValueError("Expected quadratic edges and a finite two-component line traction")
    lengths = np.linalg.norm(mesh.nodes[edges[:, 2]]-mesh.nodes[edges[:, 0]], axis=1)
    vals = lengths[:, None, None]*np.array([1/6, 2/3, 1/6])[None, :, None]*traction
    force = np.zeros(mesh.ndof)
    np.add.at(force, (2*edges[:, :, None]+[0, 1]).ravel(), vals.ravel())
    return force


def _check_rigid_constraints(mesh, fixed):
    edges = mesh.cells[:, [[0, 3], [3, 1], [1, 4], [4, 2], [2, 5], [5, 0]]].reshape(-1, 2)
    adjacency = sparse.coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                                  shape=(len(mesh.nodes), len(mesh.nodes))).tocsr()
    ncomp, labels = connected_components(adjacency, directed=False)
    for component in range(ncomp):
        selected = fixed[labels[fixed//2] == component]
        if len(selected) < 3:
            raise ValueError("Unconstrained component: rigid-body motion remains")
        xy = mesh.nodes[selected//2]
        # Center and scale rotation before rank checking to avoid unit effects.
        xy = xy-xy.mean(axis=0)
        xy /= max(float(np.max(np.abs(xy))), 1.)
        rbm = np.column_stack([selected % 2 == 0, selected % 2 == 1,
                               np.where(selected % 2 == 0, -xy[:, 1], xy[:, 0])])
        if np.linalg.matrix_rank(rbm, tol=1e-12) < 3:
            raise ValueError("Insufficient displacement constraints for rigid-body modes")


def solve_p2(mesh, moduli, material, forces=None, prescribed=None,
             body_force=None, quadrature_order=4):
    """Solve using DOF:value prescribed data and one reusable factorization.

    Prescribed values may be scalars or vectors with one value per load case.
    Inputs are independent of primary-solver boundary-condition interfaces.
    """
    started = perf_counter()
    stiffness = assemble_p2(mesh, moduli, material)
    if forces is None:
        forces = np.zeros((mesh.ndof, 1))
    forces = np.array(forces, dtype=float, copy=True)
    if forces.ndim == 1:
        forces = forces[:, None]
    if forces.ndim != 2 or forces.shape[0] != mesh.ndof or forces.shape[1] == 0 or not np.isfinite(forces).all():
        raise ValueError("Forces must be finite, with one row per DOF")
    if body_force is not None:
        forces += body_force_p2(mesh, material, body_force, quadrature_order)[:, None]
    prescribed = {} if prescribed is None else prescribed
    fixed = np.array(sorted(prescribed), dtype=int)
    if (any(not isinstance(k, (int, np.integer)) for k in prescribed)
            or np.any(fixed < 0) or np.any(fixed >= mesh.ndof)):
        raise ValueError("Invalid prescribed DOF")
    _check_rigid_constraints(mesh, fixed)
    free = np.setdiff1d(np.arange(mesh.ndof), fixed)
    u = np.zeros_like(forces)
    for dof in fixed:
        u[dof] = prescribed[dof]
    if not np.isfinite(u).all():
        raise ValueError("Prescribed displacements must be finite")
    skew = stiffness-stiffness.T
    symmetry = float(sparse.linalg.norm(skew)/max(sparse.linalg.norm(stiffness), 1.))
    if symmetry > 1e-12 or not np.isfinite(stiffness.data).all():
        raise ArithmeticError("Independent stiffness failed symmetry/finiteness check")
    diagonal_ratio = None
    if len(free):
        reduced = stiffness[free][:, free].tocsc()
        rhs = forces[free]-stiffness[free][:, fixed]@u[fixed]
        try:
            factor = splu(reduced, permc_spec="COLAMD")
            pivots = np.abs(factor.U.diagonal())
            diagonal_ratio = float(pivots.min()/pivots.max())
            # A successful LU alone does not distinguish a singular hinge from
            # roundoff-sized pivots. Flag such systems as unresolved, without
            # interpreting this ordering-dependent ratio as a condition number.
            if diagonal_ratio < 1e-14:
                raise ValueError("P2 system is singular or an unresolved near mechanism (tiny LU pivot)")
            u[free] = factor.solve(rhs)
        except RuntimeError as exc:
            raise ValueError("P2 stiffness is singular or factorization failed") from exc
    if not np.isfinite(u).all():
        raise ArithmeticError("P2 displacement contains nonfinite values")
    residual = stiffness@u-forces
    free_error = np.linalg.norm(residual[free], axis=0)/np.maximum(np.linalg.norm(forces[free], axis=0), 1.)
    if np.any(free_error > 1e-8):
        raise ArithmeticError("P2 free-DOF equilibrium residual exceeds 1e-8")
    reactions = np.zeros_like(forces)
    reactions[fixed] = residual[fixed]
    total = forces+reactions
    resultant = total.reshape(-1, 2, forces.shape[1]).sum(axis=0)
    moments = (mesh.nodes[:, 0, None]*total[1::2]-mesh.nodes[:, 1, None]*total[0::2]).sum(axis=0)
    compliance = np.einsum("ij,ij->j", forces, u)
    internal = np.einsum("ij,ij->j", u, stiffness@u)
    work = np.einsum("ij,ij->j", forces+reactions, u)
    energy_error = np.abs(work-internal)/np.maximum(np.abs(internal), 1e-12)
    gram = forces.T@u
    applied = forces.reshape(-1, 2, forces.shape[1]).sum(axis=0)
    force_scale = np.maximum(np.linalg.norm(applied, axis=0), 1.)
    characteristic_length = max(float(np.ptp(mesh.nodes, axis=0).max()), 1.)
    force_balance_error = np.linalg.norm(resultant, axis=0)/force_scale
    moment_balance_error = np.abs(moments)/np.maximum(force_scale*characteristic_length, 1.)
    if (not np.isfinite(gram).all() or not np.isfinite(internal).all()
            or np.any(energy_error > 1e-8)):
        raise ArithmeticError("P2 energy/finiteness check failed")
    if np.any(force_balance_error > 1e-7) or np.any(moment_balance_error > 1e-7):
        raise ArithmeticError("P2 reaction force or moment equilibrium check failed")
    nonzero_load = np.linalg.norm(forces, axis=0) > 0
    if np.all(u[fixed] == 0) and np.any(compliance[nonzero_load] <= 0):
        raise ArithmeticError("P2 nonzero homogeneous load has nonpositive compliance")
    diagnostics = {
        "passed": True, "solver_converged": True, "element": "P2_triangle", "triangles": len(mesh.cells),
        "nodes": len(mesh.nodes), "ndof": mesh.ndof, "stiffness_nnz": stiffness.nnz,
        "stiffness_quadrature": "three-point degree-two exact", "symmetry_error": symmetry,
        "free_residual_relative": free_error.tolist(), "energy_error_relative": energy_error.tolist(),
        "applied_resultant": forces.reshape(-1, 2, forces.shape[1]).sum(axis=0).tolist(),
        "reaction_resultant": reactions.reshape(-1, 2, forces.shape[1]).sum(axis=0).tolist(),
        "reaction_moment": (mesh.nodes[:, 0, None]*reactions[1::2]-mesh.nodes[:, 1, None]*reactions[0::2]).sum(axis=0).tolist(),
        "applied_moment": (mesh.nodes[:, 0, None]*forces[1::2]-mesh.nodes[:, 1, None]*forces[0::2]).sum(axis=0).tolist(),
        "resultant_balance": resultant.tolist(), "moment_balance": moments.tolist(),
        "resultant_balance_relative": force_balance_error.tolist(),
        "moment_balance_relative": moment_balance_error.tolist(),
        "twice_strain_energy": internal.tolist(), "prescribed_dofs": len(fixed),
        "lu_pivot_min_max_ratio": diagonal_ratio,
        "factorizations": int(bool(len(free))), "rhs": forces.shape[1] if len(free) else 0,
        "conditioning_note": "LU pivot ratio is diagnostic, not a matrix condition estimate",
    }
    return P2Result(mesh, u, forces, compliance, gram, diagnostics, perf_counter()-started)


def _boundary_patch(mesh, axis, coordinate, low, high):
    # Reconstruct exact required exterior intervals; no primary boundary tags.
    edges = mesh.boundary_edges
    points = mesh.nodes[edges]
    span = points[:, [0, 2], 1-axis]
    minimum, maximum = span.min(axis=1), span.max(axis=1)
    select = (np.all(np.isclose(points[:, :, axis], coordinate, rtol=0, atol=1e-10), axis=1)
              & (minimum >= low-1e-10) & (maximum <= high+1e-10))
    intervals = np.column_stack([minimum[select], maximum[select]])
    intervals = intervals[np.argsort(intervals[:, 0])]
    if (not len(intervals) or abs(intervals[0, 0]-low) > 1e-10
            or abs(intervals[-1, 1]-high) > 1e-10
            or not np.allclose(intervals[1:, 0], intervals[:-1, 1], rtol=0, atol=1e-10)):
        raise ValueError("Required independent support/load boundary patch is incomplete")
    return edges[select]


def bracket_verify(mesh, moduli, material, load=1000.0):
    """Verify right/down resultant-load bases on the exact bracket boundaries."""
    started = perf_counter()
    if not np.isfinite(load) or load <= 0:
        raise ValueError("Load must be positive and finite")
    p2 = p2_mesh_from_quads(mesh)
    loaded = _boundary_patch(p2, 0, 120., 10., 30.)
    supported = _boundary_patch(p2, 1, 80., 0., 40.)
    forces = np.column_stack([line_force_p2(p2, loaded, [load/20, 0]),
                              line_force_p2(p2, loaded, [0, -load/20])])
    fixed = (2*np.unique(supported)[:, None]+[0, 1]).ravel()
    result = solve_p2(p2, moduli, material, forces, {int(i): 0. for i in fixed})
    x, y = p2.nodes.T
    pad = (((x >= 0) & (x <= 40) & (y >= 70) & (y <= 80))
           | ((x >= 110) & (x <= 120) & (y >= 10) & (y <= 30)))
    nominal = result.u[:, 1].reshape(-1, 2)
    result.max_pad_displacement = float(np.linalg.norm(nominal[pad], axis=1).max())
    result.diagnostics["max_pad_displacement"] = result.max_pad_displacement
    result.diagnostics["load_patch_mean_displacement"] = (forces[:, 1]@result.u[:, 1]/load)
    result.wall_seconds = perf_counter()-started
    return result


def error_norms_p2(mesh, u, material, exact, exact_gradient, quadrature_order=4):
    """Absolute displacement L2 and strain-energy norms of a smooth fixture."""
    u = np.asarray(u, dtype=float)
    if u.shape == (mesh.ndof, 1):
        u = u[:, 0]
    if u.shape != (mesh.ndof,):
        raise ValueError("Error norms require exactly one displacement field")
    ue = u.reshape(-1, 2)[mesh.cells]
    corners = mesh.nodes[mesh.cells[:, :3]]
    law = _law(material.E, material.nu)
    l2sq, energysq = 0., 0.
    for r, s, weight in triangle_quadrature(quadrature_order):
        vals, deriv = shape_p2(r, s)
        _, gradient, det = _kinematics(corners, deriv)
        xy = (1-r-s)*corners[:, 0]+r*corners[:, 1]+s*corners[:, 2]
        diff = np.einsum("a,eai->ei", vals, ue)-exact(xy)
        gd = np.einsum("eaj,eai->eij", gradient, ue)-exact_gradient(xy)
        strain = np.column_stack([gd[:, 0, 0], gd[:, 1, 1], gd[:, 0, 1]+gd[:, 1, 0]])
        l2sq += float(np.sum(np.sum(diff*diff, axis=1)*det*weight))
        energysq += float(np.sum(np.einsum("ei,ij,ej->e", strain, law, strain)*det*weight*material.thickness))
    return {"l2": float(np.sqrt(l2sq)), "energy": float(np.sqrt(energysq))}


def manufactured_fields(material, length=100., height=40., a=.01, b=.02):
    """Declared sine solution and -div(sigma), independently differentiated."""
    kx, ky = np.pi/length, np.pi/height
    lam = material.E*material.nu/(1-material.nu**2)
    mu = material.E/(2*(1+material.nu))
    def exact(xy):
        phi = np.sin(kx*xy[:, 0])*np.sin(ky*xy[:, 1])
        return phi[:, None]*[a, b]
    def gradient(xy):
        dx = kx*np.cos(kx*xy[:, 0])*np.sin(ky*xy[:, 1])
        dy = ky*np.sin(kx*xy[:, 0])*np.cos(ky*xy[:, 1])
        return np.array([a, b])[None, :, None]*np.column_stack([dx, dy])[:, None, :]
    def force(xy):
        phi = np.sin(kx*xy[:, 0])*np.sin(ky*xy[:, 1])
        chi = np.cos(kx*xy[:, 0])*np.cos(ky*xy[:, 1])
        fx = a*((lam+2*mu)*kx*kx+mu*ky*ky)*phi-(lam+mu)*b*kx*ky*chi
        fy = b*(mu*kx*kx+(lam+2*mu)*ky*ky)*phi-(lam+mu)*a*kx*ky*chi
        return np.column_stack([fx, fy])
    return exact, gradient, force
