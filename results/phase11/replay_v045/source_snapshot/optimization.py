"""SIMP/OC baseline, pad-aware density filtering and binary geometry checks.

This is an independent implementation of the equations frozen in
research/MATHEMATICAL_FOUNDATION.md, sections 3--4. SIMP, density filtering and
optimality criteria are established methods; algorithmic attribution is to
Andreassen et al. (2011), doi:10.1007/s00158-010-0594-7. No top88 source is copied.
All volume calculations use physical densities after the solid-pad override.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

from .geometry import bracket_pads, bracket_fixed_dofs, bracket_load_basis
from . import fem


def _vector(value, size, name):
    result = np.asarray(value, dtype=float)
    if result.shape != (size,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite vector of length {size}")
    return result


class DensityFilter:
    """Area-weighted compact linear filter over active cells only.

    H maps raw values to filtered values. Preserved inputs are set to one,
    then preserved outputs are overridden to one. Consequently the derivative
    is H.T with preserved output rows and input columns removed.
    """

    def __init__(self, mesh, pads, radius):
        self.mesh = mesh
        self.pads = np.asarray(pads).copy()
        self.areas = np.asarray(mesh.areas, dtype=float)
        self.size = len(mesh.cells)
        if self.pads.dtype != bool or self.pads.shape != (self.size,):
            raise ValueError("Pads must be a Boolean vector with one value per cell")
        if not np.isfinite(radius) or radius <= 0:
            raise ValueError("Filter radius must be positive and finite")
        if not np.isfinite(self.areas).all() or np.any(self.areas <= 0):
            raise ValueError("Filter cells must have positive, finite areas")
        self.radius = float(radius)
        self.free = ~self.pads
        centers = mesh.centroids
        neighbors = cKDTree(centers).query_ball_point(centers, radius, return_sorted=True)
        lengths = np.fromiter((len(n) for n in neighbors), dtype=np.int64)
        rows = np.repeat(np.arange(self.size), lengths)
        columns = np.concatenate(neighbors)
        weights = np.maximum(0, radius-np.linalg.norm(centers[rows]-centers[columns], axis=1))
        weights *= self.areas[columns]
        totals = np.bincount(rows, weights=weights, minlength=self.size)
        if np.any(totals <= 0):
            raise ValueError("A filter row has no positive weight")
        self.H = csr_matrix((weights/totals[rows], (rows, columns)),
                            shape=(self.size, self.size))
        self.H.eliminate_zeros()

    def physical(self, raw):
        raw = _vector(raw, self.size, "Raw density").copy()
        if np.any(raw < 0) or np.any(raw > 1):
            raise ValueError("Raw density must lie in [0,1]")
        raw[self.pads] = 1.0
        physical = np.asarray(self.H @ raw)
        # A normalized nonnegative filter is a convex combination. Remove only
        # floating-point roundoff outside its mathematical [0,1] range.
        physical = np.clip(physical, 0.0, 1.0)
        physical[self.pads] = 1.0
        return physical

    def volume(self, raw, thickness):
        if not np.isfinite(thickness) or thickness <= 0:
            raise ValueError("Thickness must be positive and finite")
        return float(thickness * self.areas @ self.physical(raw))

    def chain(self, physical_gradient):
        gradient = _vector(physical_gradient, self.size, "Physical gradient").copy()
        gradient[self.pads] = 0.0
        raw_gradient = np.asarray(self.H.T @ gradient)
        raw_gradient[self.pads] = 0.0
        return raw_gradient

    def volume_gradient(self, thickness):
        if not np.isfinite(thickness) or thickness <= 0:
            raise ValueError("Thickness must be positive and finite")
        return self.chain(thickness*self.areas)

    def initialize(self, fraction=.4, xmin=.001):
        """Constant free raw values satisfying the actual filtered volume.

        The feasible endpoint of bisection is returned. A budget below the
        lower-bound filtered volume is impossible and is explicitly rejected.
        """
        if not np.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError("Volume fraction must lie in (0,1]")
        if not np.isfinite(xmin) or not 0 <= xmin <= 1:
            raise ValueError("Raw lower bound must lie in [0,1]")
        target = fraction*self.areas.sum()
        raw = np.full(self.size, xmin)
        raw[self.pads] = 1.0
        minimum = self.volume(raw, 1.0)
        if minimum > target+1e-12*self.areas.sum():
            raise ValueError("Budget is below the minimum pad-aware filtered volume")
        if fraction == 1:
            return np.ones(self.size)
        lo, hi = float(xmin), 1.0
        for _ in range(100):
            middle = .5*(lo+hi)
            raw[self.free] = middle
            if self.volume(raw, 1.0) > target:
                hi = middle
            else:
                lo = middle
            if hi-lo <= 1e-14:
                break
        raw[self.free] = lo
        return raw


def _validate_parameters(material, benchmark):
    values = [material.E, material.nu, material.thickness, benchmark.volume_fraction,
              benchmark.load, benchmark.filter_radius, benchmark.p,
              benchmark.emin_ratio, benchmark.xmin, benchmark.move]
    if not np.isfinite(values).all():
        raise ValueError("SIMP settings and material parameters must be finite")
    if material.E <= 0 or material.thickness <= 0 or not -1 < material.nu < .5:
        raise ValueError("Invalid elastic material or thickness")
    if not 0 < benchmark.volume_fraction <= 1 or benchmark.load <= 0:
        raise ValueError("Positive load and volume fraction in (0,1] required")
    if benchmark.filter_radius <= 0 or benchmark.p < 1:
        raise ValueError("Positive filter radius and SIMP exponent >=1 required")
    if not 0 < benchmark.emin_ratio < 1 or not 0 < benchmark.xmin < 1:
        raise ValueError("Emin ratio and raw lower bound must lie in (0,1)")
    if not 0 < benchmark.move <= 1:
        raise ValueError("OC move limit must lie in (0,1]")
    for name in ("max_evaluations", "min_evaluations"):
        value = getattr(benchmark, name)
        if isinstance(value, bool) or int(value) != value or value < 1:
            raise ValueError(f"{name} must be a positive integer")


def compliance_sensitivity(mesh, physical, displacement, material, benchmark,
                           unit_stiffness=None):
    """Derivative dC/d(rho), for the single nominal load solve.

    The density-filter transpose and pad override are applied separately by
    DensityFilter.chain. Stiffness includes thickness and unit Young modulus.
    """
    rho = _vector(physical, len(mesh.cells), "Physical density")
    if np.any(rho < 0) or np.any(rho > 1):
        raise ValueError("Physical density must lie in [0,1]")
    u = np.asarray(displacement, dtype=float)
    if u.shape == (mesh.ndof, 1):
        u = u[:, 0]
    u = _vector(u, mesh.ndof, "Displacement")
    if unit_stiffness is None:
        unit_stiffness = fem.unit_stiffness(mesh, material.nu, material.thickness)
    ke = np.asarray(unit_stiffness)
    if ke.shape != (len(mesh.cells), 8, 8) or not np.isfinite(ke).all():
        raise ValueError("Unit stiffness must have shape (ne,8,8) and finite values")
    edofs = (2*mesh.cells[:, :, None]+np.array([0, 1])).reshape(-1, 8)
    ue = u[edofs]
    energy = np.einsum("ei,eij,ej->e", ue, ke, ue, optimize=True)
    return -benchmark.p*rho**(benchmark.p-1)*material.E*(1-benchmark.emin_ratio)*energy


def oc_update(raw, compliance_gradient, density_filter, thickness,
              fraction=.4, xmin=.001, move=.2):
    """Return the feasible OC bisection endpoint; never evaluate an FE solve.

    Every multiplier test uses filtered physical volume. A zero gradient is
    permitted. Substantial positive sensitivities or impossible moves fail.
    """
    size = density_filter.size
    x = _vector(raw, size, "Raw density").copy()
    g = _vector(compliance_gradient, size, "Compliance gradient")
    free = density_filter.free
    if (not np.isfinite([fraction, xmin, move]).all() or not 0 < fraction <= 1
            or not 0 <= xmin < 1 or not 0 < move <= 1):
        raise ValueError("Invalid OC fraction, lower bound, or move limit")
    if np.any(x[free] < xmin) or np.any(x > 1) or np.any(x[density_filter.pads] != 1):
        raise ValueError("OC raw field violates bounds or fixed pads")
    h = density_filter.volume_gradient(thickness)
    if np.any(h[free] <= 0):
        raise ValueError("Every free volume sensitivity must be positive")
    scale = max(1.0, float(np.max(np.abs(g[free]), initial=0)))
    if np.any(g[free] > 1e-10*scale):
        raise ValueError("Compliance has a substantially positive sensitivity")
    target = fraction*thickness*density_filter.areas.sum()
    tolerance = 1e-8*thickness*density_filter.areas.sum()
    lower = np.maximum(xmin, x-move)
    upper = np.minimum(1.0, x+move)
    lower[density_filter.pads] = upper[density_filter.pads] = 1.0
    if density_filter.volume(lower, thickness) > target+1e-12*max(target, 1):
        raise ValueError("No volume-feasible OC move-limit candidate exists")
    if density_filter.volume(upper, thickness) <= target:
        return upper

    def candidate(multiplier):
        candidate_x = x.copy()
        candidate_x[free] = np.clip(x[free]*np.sqrt(np.maximum(0, -g[free]) /
                                                 (multiplier*h[free])),
                                    lower[free], upper[free])
        return candidate_x

    lo, hi = 0.0, 1.0
    feasible = candidate(hi)
    while density_filter.volume(feasible, thickness) > target:
        hi *= 2
        if not np.isfinite(hi):
            raise ValueError("Could not bracket a feasible OC multiplier")
        feasible = candidate(hi)
    for _ in range(100):
        middle = .5*(lo+hi)
        trial = candidate(middle)
        trial_volume = density_filter.volume(trial, thickness)
        if trial_volume > target:
            lo = middle
        else:
            hi, feasible = middle, trial
        gap = target-density_filter.volume(feasible, thickness)
        if gap <= tolerance or (hi-lo)/max(hi+lo, np.finfo(float).tiny) <= 1e-8:
            break
    return feasible


@dataclass
class SimpResult:
    raw: np.ndarray
    physical: np.ndarray
    history: list
    status: str
    evaluations: int
    solve_seconds: float
    best_raw: np.ndarray
    best_physical: np.ndarray
    best_evaluation: int


def solve_simp(mesh, material, benchmark, initial=None, callback=None):
    """Run the deterministic nominal-load SIMP/OC baseline.

    The initial state consumes one evaluation. A stop decision precedes every
    update, so raw/physical always identify the last *evaluated* state. The
    lowest-compliance feasible evaluated incumbent is also returned. Solver
    failures raise with optimization history and candidate arrays attached.
    """
    _validate_parameters(material, benchmark)
    pads = bracket_pads(mesh)
    filt = DensityFilter(mesh, pads, benchmark.filter_radius)
    if initial is None:
        raw = filt.initialize(benchmark.volume_fraction, benchmark.xmin)
    else:
        raw = _vector(initial, len(mesh.cells), "Initial density").copy()
        if (np.any(raw[~pads] < benchmark.xmin) or np.any(raw > 1)
                or np.any(raw[pads] != 1)):
            raise ValueError("Initial density violates bounds or preserved pads")
    full_volume = float(material.thickness*mesh.areas.sum())
    budget = benchmark.volume_fraction*full_volume
    if filt.volume(raw, material.thickness) > budget+1e-6*full_volume:
        raise ValueError("Initial physical volume exceeds the budget")
    force = benchmark.load*bracket_load_basis(mesh)[:, 1]
    fixed = bracket_fixed_dofs(mesh)
    unit_ke = fem.unit_stiffness(mesh, material.nu, material.thickness)
    history = []
    previous_raw = None
    best_raw, best_physical = None, None
    best_compliance, best_evaluation = np.inf, 0
    solve_seconds = 0.0
    status = "budget_exhausted"
    for evaluation in range(1, benchmark.max_evaluations+1):
        physical = filt.physical(raw)
        volume = float(material.thickness*mesh.areas @ physical)
        if volume > budget+1e-6*full_volume:
            raise ValueError("An evaluated candidate exceeds physical volume tolerance")
        moduli = material.E*(benchmark.emin_ratio+
                             (1-benchmark.emin_ratio)*physical**benchmark.p)
        try:
            solved = fem.solve(mesh, moduli, force, fixed, material.nu, material.thickness)
            objective = float(solved.compliance[0])
            gradient = filt.chain(compliance_sensitivity(mesh, physical, solved.u,
                                                        material, benchmark, unit_ke))
            if not np.isfinite(gradient).all() or not np.isfinite(objective) or objective <= 0:
                raise ValueError("Nonfinite sensitivity or nonpositive compliance")
        except Exception as exc:
            exc.optimization_history = history
            exc.last_evaluated_raw = None if previous_raw is None else previous_raw.copy()
            exc.failed_raw = raw.copy()
            raise
        solve_seconds += solved.wall_seconds
        max_change = None if previous_raw is None else float(np.max(np.abs(raw-previous_raw)))
        compliance_window = [row["compliance"] for row in history[-10:]]+[objective]
        spread = ((max(compliance_window)-min(compliance_window))/max(abs(objective), 1e-12)
                  if len(compliance_window) == 11 else None)
        row = {"evaluation": evaluation, "compliance": objective,
               "volume": volume, "volume_fraction": volume/full_volume,
               "volume_violation": max(0.0, volume-budget),
               "max_change": max_change, "relative_compliance_spread": spread,
               "gradient_min": float(np.min(gradient[~pads])),
               "gradient_max": float(np.max(gradient[~pads])),
               "solve_seconds": float(solved.wall_seconds),
               "cumulative_solve_seconds": float(solve_seconds),
               "diagnostics": solved.diagnostics}
        history.append(row)
        if objective < best_compliance:
            best_compliance, best_evaluation = objective, evaluation
            best_raw, best_physical = raw.copy(), physical.copy()
        if callback is not None:
            callback(row)
        changes = [r["max_change"] for r in history[-10:]]
        converged = (evaluation >= benchmark.min_evaluations and len(history) >= 11
                     and all(change is not None and change <= .01 for change in changes)
                     and spread <= .001)
        if converged:
            status = "converged"
            break
        if evaluation == benchmark.max_evaluations:
            break
        previous_raw = raw.copy()
        try:
            raw = oc_update(raw, gradient, filt, material.thickness,
                            benchmark.volume_fraction, benchmark.xmin, benchmark.move)
        except Exception as exc:
            exc.optimization_history = history
            exc.last_evaluated_raw = previous_raw.copy()
            raise
    return SimpResult(raw.copy(), physical.copy(), history, status, len(history),
                      solve_seconds, best_raw, best_physical, best_evaluation)


def _binary_count(mesh, fraction, pads):
    if not np.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError("Binary fraction must lie in (0,1]")
    areas = mesh.areas
    if not np.allclose(areas, areas[0], rtol=0, atol=1e-10):
        raise ValueError("Rank conversion requires equal-area cells")
    desired = fraction*len(mesh.cells)
    count = round(desired)
    if not np.isclose(desired, count, rtol=0, atol=1e-9):
        raise ValueError("The requested binary volume is not an integer cell count")
    if count < np.count_nonzero(pads):
        raise ValueError("Binary budget is smaller than preserved pads")
    return count-int(np.count_nonzero(pads))


def _rank_binary(mesh, score, fraction):
    pads = bracket_pads(mesh)
    fill = _binary_count(mesh, fraction, pads)
    free = np.flatnonzero(~pads)
    xy = mesh.cell_indices
    order = np.lexsort((xy[free, 0], xy[free, 1], score[free]))
    binary = pads.copy()
    binary[free[order[:fill]]] = True
    return binary


def simple_l_design(mesh, fraction=.4):
    """B1: nearest centroids to (20,75)--(20,20)--(115,20), then y,x.

    At h=2 and fraction=.4 this selects 150 pad cells and 490 free cells.
    Frozen geometry is refined by transferring this field, never reranking.
    """
    centers = mesh.centroids
    starts = np.array([[20., 75.], [20., 20.]])
    ends = np.array([[20., 20.], [115., 20.]])
    distances = []
    for start, end in zip(starts, ends):
        tangent = end-start
        parameter = np.clip((centers-start) @ tangent/(tangent @ tangent), 0, 1)
        nearest = start+parameter[:, None]*tangent
        distances.append(np.linalg.norm(centers-nearest, axis=1))
    return _rank_binary(mesh, np.minimum(*distances), fraction)


def binary_by_rank(mesh, physical, fraction=.4):
    """Retain fixed pads and highest free physical densities, ties by y,x."""
    rho = _vector(physical, len(mesh.cells), "Physical density")
    if np.any(rho < 0) or np.any(rho > 1):
        raise ValueError("Physical density must lie in [0,1]")
    return _rank_binary(mesh, -rho, fraction)


def geometry_status(mesh, binary):
    """Check pad presence and edge connectivity, rejecting vertex-only joins.

    No disconnected component is removed or restrained. Offending area is the
    total solid area outside the first support-connected component.
    """
    binary = np.asarray(binary)
    if binary.shape != (len(mesh.cells),) or not np.isin(binary, [0, 1]).all():
        raise ValueError("Geometry requires one binary value per element")
    solid = binary.astype(bool)
    pads = bracket_pads(mesh)
    centers = mesh.centroids
    support_pad = pads & (centers[:, 1] >= 70)
    load_pad = pads & (centers[:, 0] >= 110)
    support_present = bool(np.isclose(mesh.areas[support_pad].sum(), 400., rtol=0,
                                      atol=1e-8) and np.all(solid[support_pad]))
    load_present = bool(np.isclose(mesh.areas[load_pad].sum(), 200., rtol=0,
                                   atol=1e-8) and np.all(solid[load_pad]))
    active = np.flatnonzero(solid)
    neighbors = {int(cell): [] for cell in active}
    owner = {}
    for cell in active:
        vertices = mesh.cells[cell]
        for j in range(4):
            edge = tuple(sorted((int(vertices[j]), int(vertices[(j+1) % 4]))))
            if edge in owner:
                other = owner[edge]
                neighbors[int(cell)].append(other)
                neighbors[other].append(int(cell))
            else:
                owner[edge] = int(cell)
    labels = np.full(len(mesh.cells), -1, dtype=int)
    components = []
    for seed in active:
        if labels[seed] >= 0:
            continue
        component = []
        stack = [int(seed)]
        label = len(components)
        labels[seed] = label
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in neighbors[current]:
                if labels[neighbor] < 0:
                    labels[neighbor] = label
                    stack.append(neighbor)
        components.append(component)
    support_cells = np.flatnonzero(solid & support_pad)
    support_label = int(labels[support_cells[0]]) if len(support_cells) else -1
    offending = solid & (labels != support_label) if support_label >= 0 else solid
    pads_connected = bool(support_present and load_present and
                          np.all(labels[pads] == support_label))
    connected = bool(len(components) == 1 and pads_connected)
    return {"valid": connected, "connected": connected,
            "pads_preserved": support_present and load_present,
            "pads_present": support_present and load_present,
            "support_pad_present": support_present, "load_pad_present": load_present,
            "pads_connected": pads_connected, "component_count": len(components),
            "solid_cells": int(solid.sum()), "solid_area": float(mesh.areas @ solid),
            "offending_area": float(mesh.areas @ offending),
            "component_areas": [float(mesh.areas[c].sum()) for c in components]}


def sensitivity_audit(mesh, material, benchmark, seed=20261010, signed=False):
    """Ten reproducible central-difference sweeps, including pad/cutout edges.

    Returns numeric arrays for archival as NPZ. Directions have zero pad entries
    and unit Euclidean norm. No perturbed evaluation is clipped or normalized.
    The first three directions are localized near support pad, load pad, and
    re-entrant corner respectively; seven are globally distributed directions.
    Nonnegative entries avoid cancellation of negative compliance sensitivities
    and keep the directional derivative well scaled relative to FE solve error.
    Set signed=True for seven additional global sign-changing directions; their
    cancellation may obscure a truncation-convergence region at these steps.
    This diagnostic is separate from and charged outside the optimizer budget.
    """
    _validate_parameters(material, benchmark)
    filt = DensityFilter(mesh, bracket_pads(mesh), benchmark.filter_radius)
    rng = np.random.default_rng(seed)
    raw = .35+.3*rng.random(len(mesh.cells))
    raw[filt.pads] = 1.0
    directions = rng.normal(size=(10, len(raw)))
    if not signed:
        directions = np.abs(directions)
    for index, center in enumerate(((20., 68.), (108., 20.), (40., 40.))):
        distances = np.linalg.norm(mesh.centroids-np.array(center), axis=1)
        directions[index, distances > max(12.0, 2*mesh.h)] = 0.0
        directions[index] = np.abs(directions[index])
    directions[:, filt.pads] = 0.0
    norms = np.linalg.norm(directions, axis=1)
    if np.any(norms == 0):
        raise ValueError("Mesh does not resolve all three localized audit directions")
    directions /= norms[:, None]
    steps = np.logspace(-2, -6, 5)
    force = benchmark.load*bracket_load_basis(mesh)[:, 1]
    fixed = bracket_fixed_dofs(mesh)
    solve_seconds = 0.0

    def evaluate(x):
        nonlocal solve_seconds
        rho = filt.physical(x)
        modulus = material.E*(benchmark.emin_ratio+(1-benchmark.emin_ratio)*rho**benchmark.p)
        result = fem.solve(mesh, modulus, force, fixed, material.nu, material.thickness)
        solve_seconds += result.wall_seconds
        return result, rho

    baseline, physical = evaluate(raw)
    gradient = filt.chain(compliance_sensitivity(mesh, physical, baseline.u, material, benchmark))
    volume_gradient = filt.volume_gradient(material.thickness)
    exact_compliance = directions @ gradient
    exact_volume = directions @ volume_gradient
    compliance_fd = np.zeros((len(directions), len(steps)))
    volume_fd = np.zeros_like(compliance_fd)
    for index, direction in enumerate(directions):
        for step_index, step in enumerate(steps):
            plus, minus = raw+step*direction, raw-step*direction
            plus_result, _ = evaluate(plus)
            minus_result, _ = evaluate(minus)
            compliance_fd[index, step_index] = (plus_result.compliance[0]-minus_result.compliance[0])/(2*step)
            volume_fd[index, step_index] = (filt.volume(plus, material.thickness)-
                                           filt.volume(minus, material.thickness))/(2*step)
    compliance_absolute = np.abs(compliance_fd-exact_compliance[:, None])
    volume_absolute = np.abs(volume_fd-exact_volume[:, None])
    compliance_scale = np.maximum(np.maximum(np.abs(compliance_fd),
                                             np.abs(exact_compliance[:, None])), 1e-10)
    volume_scale = np.maximum(np.maximum(np.abs(volume_fd), np.abs(exact_volume[:, None])), 1e-10)
    compliance_error = compliance_absolute/compliance_scale
    volume_error = volume_absolute/volume_scale
    passed_steps = compliance_error < 1e-4
    adjacent = np.any(passed_steps[:, :-1] & passed_steps[:, 1:], axis=1)
    return {"seed": np.array(seed), "signed_global_directions": np.array(signed),
            "raw": raw, "physical": physical,
            "directions": directions, "steps": steps,
            "compliance_gradient": gradient, "volume_gradient": volume_gradient,
            "compliance_directional_derivatives": exact_compliance,
            "volume_directional_derivatives": exact_volume,
            "compliance_fd": compliance_fd, "volume_fd": volume_fd,
            "compliance_absolute_errors": compliance_absolute,
            "compliance_relative_errors": compliance_error,
            "volume_absolute_errors": volume_absolute,
            "volume_relative_errors": volume_error,
            "compliance_two_adjacent_pass": adjacent,
            "factorizations": np.array(1+2*len(directions)*len(steps)),
            "solve_seconds": np.array(solve_seconds)}
