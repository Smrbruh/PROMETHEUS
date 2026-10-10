"""Deterministic conforming square meshes and geometric boundary integration.

Coordinates and cell size are in mm. Connectivity is counterclockwise
bottom-left, bottom-right, top-right, top-left; cells are ordered by y then x.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Mesh:
    nodes: np.ndarray
    cells: np.ndarray
    cell_indices: np.ndarray
    h: float

    def __post_init__(self):
        nodes = np.array(self.nodes, dtype=float, copy=True)
        original_cells = np.asarray(self.cells)
        original_indices = np.asarray(self.cell_indices)
        if not np.issubdtype(original_cells.dtype, np.integer):
            raise ValueError("Cell connectivity must contain integer node IDs")
        if not np.issubdtype(original_indices.dtype, np.integer):
            raise ValueError("Cell indices must be integers")
        cells = original_cells.astype(np.int64, copy=True)
        indices = original_indices.astype(np.int64, copy=True)
        if nodes.ndim != 2 or nodes.shape[1] != 2 or not np.isfinite(nodes).all():
            raise ValueError("Nodes must be a finite (n,2) array")
        if cells.ndim != 2 or cells.shape[1] != 4 or len(cells) == 0:
            raise ValueError("Mesh must contain four-node cells")
        if indices.shape != (len(cells), 2):
            raise ValueError("Cell indices must have shape (ne,2)")
        if np.any(cells < 0) or np.any(cells >= len(nodes)):
            raise ValueError("Connectivity contains an out-of-range node")
        if np.any(np.diff(np.sort(cells, axis=1), axis=1) == 0):
            raise ValueError("Element contains repeated nodes")
        if len(np.unique(cells, axis=0)) != len(cells):
            raise ValueError("Mesh contains duplicate cells")
        if not np.isfinite(self.h) or self.h <= 0:
            raise ValueError("Mesh spacing must be positive and finite")
        for name, array in (("nodes", nodes), ("cells", cells), ("cell_indices", indices)):
            array.setflags(write=False)
            object.__setattr__(self, name, array)

    @property
    def centroids(self):
        return self.nodes[self.cells].mean(axis=1)

    @property
    def areas(self):
        points = self.nodes[self.cells]
        following = np.roll(points, -1, axis=1)
        return .5 * np.sum(points[:, :, 0] * following[:, :, 1]
                           - points[:, :, 1] * following[:, :, 0], axis=1)

    @property
    def ndof(self):
        return 2 * len(self.nodes)


def make_rect_mesh(width, height, h):
    """Create a regular mesh; the spacing must divide both dimensions."""
    if not np.isfinite([width, height, h]).all() or min(width, height, h) <= 0:
        raise ValueError("Dimensions and spacing must be positive and finite")
    nx, ny = round(width / h), round(height / h)
    if nx < 1 or ny < 1 or not np.allclose([nx*h, ny*h], [width, height], rtol=0, atol=1e-10):
        raise ValueError("Spacing must divide both rectangle dimensions")
    xx, yy = np.meshgrid(np.arange(nx+1)*h, np.arange(ny+1)*h)
    nodes = np.column_stack([xx.ravel(), yy.ravel()])
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny))
    bl = (iy*(nx+1)+ix).ravel()
    cells = np.column_stack([bl, bl+1, bl+nx+2, bl+nx+1])
    return Mesh(nodes, cells, np.column_stack([ix.ravel(), iy.ravel()]), float(h))


def subset_mesh(mesh, mask):
    """Retain selected cells, preserving their order, and remove unused nodes."""
    mask = np.asarray(mask)
    if mask.dtype != bool or mask.shape != (len(mesh.cells),):
        raise ValueError("Subset mask must be a Boolean per cell")
    if not mask.any():
        raise ValueError("An empty material mesh cannot be analyzed")
    retained = mesh.cells[mask]
    used, inverse = np.unique(retained, return_inverse=True)
    return Mesh(mesh.nodes[used], inverse.reshape(-1, 4), mesh.cell_indices[mask], mesh.h)


def make_l_mesh(h=2):
    """PROMETHEUS-01 domain with its cutout absent (never ersatz-filled)."""
    if not np.isfinite(h) or h <= 0 or not np.isclose(40/h, round(40/h), rtol=0, atol=1e-10):
        raise ValueError("Spacing must divide the 40 mm legs")
    mesh = make_rect_mesh(120, 80, h)
    centers = mesh.centroids
    return subset_mesh(mesh, (centers[:, 0] < 40) | (centers[:, 1] < 40))


def bracket_pads(mesh):
    centers = mesh.centroids
    x, y = centers.T
    tol = 1e-10
    support = (x >= -tol) & (x <= 40+tol) & (y >= 70-tol) & (y <= 80+tol)
    load = (x >= 110-tol) & (x <= 120+tol) & (y >= 10-tol) & (y <= 30+tol)
    return support | load


def boundary_edges(mesh):
    """Return oriented edges occurring exactly once in the cell mesh."""
    edges = mesh.cells[:, [[0, 1], [1, 2], [2, 3], [3, 0]]].reshape(-1, 2)
    _, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0,
                                   return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("Non-manifold mesh edge")
    return edges[counts[inverse] == 1]


def _patch_edges(mesh, axis, coordinate, low, high):
    edges = boundary_edges(mesh)
    xy = mesh.nodes[edges]
    along = 1-axis
    on_line = np.all(np.isclose(xy[:, :, axis], coordinate, rtol=0, atol=1e-10), axis=1)
    lo = np.min(xy[:, :, along], axis=1)
    hi = np.max(xy[:, :, along], axis=1)
    selected = on_line & (lo >= low-1e-10) & (hi <= high+1e-10)
    patch = edges[selected]
    spans = np.column_stack([lo[selected], hi[selected]])
    spans = spans[np.argsort(spans[:, 0])]
    if (len(spans) == 0 or abs(spans[0, 0]-low) > 1e-10
            or abs(spans[-1, 1]-high) > 1e-10
            or not np.allclose(spans[1:, 0], spans[:-1, 1], rtol=0, atol=1e-10)):
        raise ValueError("Required boundary patch is missing, nonconforming, or incomplete")
    return patch


def bracket_fixed_dofs(mesh):
    edges = _patch_edges(mesh, axis=1, coordinate=80, low=0, high=40)
    nodes = np.unique(edges)
    return (2*nodes[:, None] + np.array([0, 1])).ravel()


def bracket_load_basis(mesh):
    """Unit-resultant RIGHT and DOWN loads on the full 20 mm load edge.

    Integrates linear edge shape functions, including endpoint half weights.
    A line load has already absorbed thickness, so no thickness enters here.
    """
    edges = _patch_edges(mesh, axis=0, coordinate=120, low=10, high=30)
    lengths = np.linalg.norm(np.diff(mesh.nodes[edges], axis=1)[:, 0], axis=1)
    weights = np.zeros(len(mesh.nodes))
    np.add.at(weights, edges.ravel(), np.repeat(lengths/40, 2))
    basis = np.zeros((mesh.ndof, 2))
    basis[0::2, 0] = weights
    basis[1::2, 1] = -weights
    return basis
