"""Archive measured analytical and manufactured-solution validation evidence.

From the repository root:
    python scripts/validation_evidence.py --output results/validation.json

Requires the same NumPy, SciPy and pytest environment as the automated tests.
The primary Q4 MMS fixture is loaded from its tested helper with runpy; no
pytest tests are executed implicitly. The P2 analytical fields and all P2
assembly, forcing and integration use the independent verifier implementation.
An existing output file is preserved; choose another --output for a replay.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import runpy
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import scipy
from prometheus import fem
from prometheus.geometry import make_rect_mesh
from prometheus.models import Material
from prometheus.verification import (
    p2_mesh_from_quads, line_force_p2, solve_p2, manufactured_fields,
    error_norms_p2, shape_p2,
)


def analytic_strip(material):
    """1 MPa tension: free lateral contraction, one transverse pin at origin."""
    width, height, h, traction = 100., 10., 2., 1.
    mesh = make_rect_mesh(width, height, h)
    right = np.flatnonzero(mesh.nodes[:, 0] == width)
    force = np.zeros(mesh.ndof)
    force[2*right] = traction*material.thickness*h
    force[2*right[[0, -1]]] *= .5
    left = np.flatnonzero(mesh.nodes[:, 0] == 0)
    fixed = np.r_[2*left, 1]
    q4 = fem.solve(mesh, material.E, force, fixed, material.nu, material.thickness)
    q4_exact = (mesh.nodes*np.array([traction/material.E, -material.nu*traction/material.E])).ravel()
    q4_stress = []
    for cell, dofs in zip(mesh.cells, fem.element_dofs(mesh)):
        for _, bmat, _, _ in fem.quadrature_data(mesh.nodes[cell]):
            q4_stress.append(fem.plane_stress(material.E, material.nu)@bmat@q4.u[dofs, 0])
    q4_stress = np.asarray(q4_stress)

    p2mesh = p2_mesh_from_quads(mesh)
    right_edges = p2mesh.boundary_edges[
        np.all(p2mesh.nodes[p2mesh.boundary_edges, 0] == width, axis=1)]
    p2force = line_force_p2(p2mesh, right_edges, [traction*material.thickness, 0.])
    left = np.flatnonzero(p2mesh.nodes[:, 0] == 0)
    p2fixed = {int(2*i): 0. for i in left}
    origin = int(np.flatnonzero(np.all(p2mesh.nodes == [0, 0], axis=1))[0])
    p2fixed[2*origin+1] = 0.
    p2 = solve_p2(p2mesh, material.E, material, p2force, p2fixed)
    p2_exact = (p2mesh.nodes*np.array([traction/material.E, -material.nu*traction/material.E])).ravel()
    # Recover engineering strain at triangle centroids, using the P2 reference
    # derivatives and independently formed physical affine Jacobians.
    corners = p2mesh.nodes[p2mesh.cells[:, :3]]
    jac = np.stack([corners[:, 1]-corners[:, 0], corners[:, 2]-corners[:, 0]], axis=2)
    gradient = np.einsum("ai,eij->eaj", shape_p2(1/3, 1/3)[1], np.linalg.inv(jac))
    nodal = p2.u[:, 0].reshape(-1, 2)[p2mesh.cells]
    displacement_gradient = np.einsum("eaj,eai->eij", gradient, nodal)
    strain = np.column_stack([displacement_gradient[:, 0, 0], displacement_gradient[:, 1, 1],
                               displacement_gradient[:, 0, 1]+displacement_gradient[:, 1, 0]])
    lam = material.E*material.nu/(1-material.nu**2)
    mu = material.E/(2*(1+material.nu))
    p2_stress = np.column_stack([(lam+2*mu)*strain[:, 0]+lam*strain[:, 1],
                                 lam*strain[:, 0]+(lam+2*mu)*strain[:, 1], mu*strain[:, 2]])
    resultant = traction*height*material.thickness
    expected_compliance = resultant*traction*width/material.E
    results = {}
    for name, result, coordinates, exact, stress in (
            ("Q4", q4, mesh.nodes, q4_exact, q4_stress),
            ("P2", p2, p2mesh.nodes, p2_exact, p2_stress)):
        relative_u = float(np.linalg.norm(result.u[:, 0]-exact)/np.linalg.norm(exact))
        relative_c = abs(float(result.compliance[0])-expected_compliance)/expected_compliance
        stress_error = float(np.max(np.abs(stress-[traction, 0., 0.])))
        results[name] = {
            "nodes": len(coordinates), "ndof": 2*len(coordinates),
            "compliance_n_mm": float(result.compliance[0]),
            "compliance_relative_error": relative_c,
            "displacement_relative_l2_nodal_error": relative_u,
            "displacement_max_absolute_error_mm": float(np.max(np.abs(result.u[:, 0]-exact))),
            "stress_max_absolute_component_error_mpa": stress_error,
            "reaction_resultant_n": (q4.reactions.reshape(-1, 2).sum(axis=0).tolist()
                                        if name == "Q4" else np.asarray(p2.diagnostics["reaction_resultant"])[:, 0].tolist()),
            "solver_diagnostics": result.diagnostics,
            "factorizations": 1, "rhs_solves": 1,
            "passed": bool(relative_u < 1e-8 and relative_c < 1e-8 and stress_error/traction < 1e-8),
        }
    return {
        "geometry_mm": {"length": width, "height": height, "h": h},
        "traction_mpa": traction, "resultant_n": [resultant, 0.],
        "constraints": "ux=0 on x=0; uy=0 only at (0,0); other lateral boundaries traction free",
        "analytical_displacement": "ux=x/E; uy=-nu*y/E for traction=1 MPa",
        "analytical_stress_mpa": [traction, 0., 0.],
        "analytical_compliance_n_mm": expected_compliance,
        "acceptance_relative_tolerance": 1e-8,
        "results": results,
    }


def manufactured(material):
    q4_helper = runpy.run_path(str(ROOT/"tests"/"test_fem.py"))["manufactured_errors"]
    q4_rows, p2_rows = [], []
    exact, gradient, body_force = manufactured_fields(material)
    for h in (5., 2.5, 1.25):
        q4_4 = q4_helper(h, order=4)
        q4_6 = q4_helper(h, order=6)
        q4_quadrature_change = {key: abs(q4_4[key]-q4_6[key])/q4_6[key]
                                for key in ("l2_error", "energy_error")}
        q4_rows.append({"h_mm": h, "order4": q4_4, "order6": q4_6,
                        "quadrature_relative_change": q4_quadrature_change,
                        "quadrature_passed": max(q4_quadrature_change.values()) < 1e-6,
                        "factorizations": 2, "rhs_solves": 2})
        p2mesh = p2_mesh_from_quads(make_rect_mesh(100., 40., h))
        boundary = np.unique(p2mesh.boundary_edges)
        fixed = {int(d): 0. for d in (2*boundary[:, None]+[0, 1]).ravel()}
        solutions = {order: solve_p2(p2mesh, material.E, material, prescribed=fixed,
                                      body_force=body_force, quadrature_order=order) for order in (4, 6)}
        norms4 = error_norms_p2(p2mesh, solutions[4].u, material, exact, gradient, 4)
        norms6 = error_norms_p2(p2mesh, solutions[6].u, material, exact, gradient, 6)
        qchange = {key: abs(norms4[key]-norms6[key])/norms6[key] for key in ("l2", "energy")}
        displacement_change = float(np.linalg.norm(solutions[4].u-solutions[6].u)/np.linalg.norm(solutions[6].u))
        p2_rows.append({"h_mm": h, "triangles": len(p2mesh.cells), "ndof": p2mesh.ndof,
                        "order4": norms4, "order6": norms6,
                        "quadrature_relative_change": qchange,
                        "quadrature_displacement_relative_change": displacement_change,
                        "quadrature_passed": max(qchange.values()) < 1e-4 and displacement_change < 1e-7,
                        "solver_order4": solutions[4].diagnostics,
                        "solver_order6": solutions[6].diagnostics,
                        "factorizations": 2, "rhs_solves": 2})
    q4_errors = np.array([[r["order6"]["l2_error"], r["order6"]["energy_error"]] for r in q4_rows])
    p2_errors = np.array([[r["order6"]["l2"], r["order6"]["energy"]] for r in p2_rows])
    q4_rates = np.log2(q4_errors[:-1]/q4_errors[1:])
    p2_rates = np.log2(p2_errors[:-1]/p2_errors[1:])
    return {
        "fixture": "research/VALIDATION_PROTOCOL.md, Gate V1 sine manufactured solution",
        "geometry_mm": {"length": 100., "height": 40.},
        "amplitudes_mm": {"a": .01, "b": .02},
        "displacement": "(a,b)*sin(pi*x/100)*sin(pi*y/40)",
        "body_force": "-div(sigma(exact_u)), including mixed derivatives; physical N/mm^3 multiplied by thickness in weak form",
        "constraints": "Exact zero displacement on every exterior boundary node",
        "error_norm_definitions": {
            "l2": "sqrt(integral_Omega ||u_h-u_exact||^2 dA)",
            "energy": "sqrt(integral_Omega thickness*(strain_error)^T D (strain_error) dA)"},
        "quadrature_note": "Q4 stiffness always2x2; P2 stiffness degree-two three-point rule. Orders4 and6 apply to body force and error integration; P2 uses a Duffy tensor Gauss rule.",
        "rates_definition": "log2(error_at_h/error_at_h_over_2), using order6 values",
        "Q4": {"rows": q4_rows, "l2_rates": q4_rates[:, 0].tolist(),
                "energy_rates": q4_rates[:, 1].tolist(), "expected_l2_rate": 2., "expected_energy_rate": 1.,
                "rate_tolerances": {"l2_absolute": .08, "energy_absolute": .04},
                "passed": bool(np.all(np.abs(q4_rates-[2., 1.]) < [.08, .04]) and all(r["quadrature_passed"] for r in q4_rows))},
        "P2": {"rows": p2_rows, "l2_rates": p2_rates[:, 0].tolist(),
                "energy_rates": p2_rates[:, 1].tolist(), "expected_l2_rate": 3., "expected_energy_rate": 2.,
                "rate_tolerances": {"l2_absolute": .15, "energy_absolute": .1},
                "passed": bool(np.all(np.abs(p2_rates-[3., 2.]) < [.15, .1]) and all(r["quadrature_passed"] for r in p2_rows))},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results"/"validation.json")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve the existing evidence; use another --output path for a replay")
    started = perf_counter()
    material = Material()
    mms = manufactured(material)
    strip = analytic_strip(material)
    source_files = [Path(__file__).resolve(), ROOT/"tests"/"test_fem.py"]
    source_files += sorted((ROOT/"prometheus").glob("*.py"))
    source_files += [ROOT/"research"/"VALIDATION_PROTOCOL.md", ROOT/"research"/"MATHEMATICAL_FOUNDATION.md"]
    data = {
        "schema_version": 1, "research_date": "2026-10-10",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv, "material": asdict(material),
        "units": {"coordinates": "mm", "displacement": "mm", "force": "N", "modulus": "MPa",
                  "compliance": "N mm", "l2_displacement_error": "mm^2", "energy_error": "sqrt(N mm)"},
        "environment": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                        "platform": platform.platform(),
                        "threads": {k: os.environ[k] for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}},
        "source_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        "manufactured_solution": mms, "analytical_tensile_strip": strip,
        "all_checks_passed": bool(mms["Q4"]["passed"] and mms["P2"]["passed"] and all(r["passed"] for r in strip["results"].values())),
        "measured_cost": {"factorizations": 14, "rhs_solves": 14, "wall_seconds": perf_counter()-started},
        "scope": "Code verification on smooth continuum fixtures. Not physical validation, bracket solution convergence, stress safety, or a certification result.",
        "independence": "P2 shape functions, Lamé constitutive construction, stiffness/body/edge-load assembly and boundary constraints are independent from Q4. Both use NumPy/SciPy, the prescribed analytical fixture, and regular rectangular geometry.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), "all_checks_passed": data["all_checks_passed"],
                      "q4_l2_rates": mms["Q4"]["l2_rates"], "q4_energy_rates": mms["Q4"]["energy_rates"],
                      "p2_l2_rates": mms["P2"]["l2_rates"], "p2_energy_rates": mms["P2"]["energy_rates"],
                      "wall_seconds": data["measured_cost"]["wall_seconds"]}, indent=2))
    if not data["all_checks_passed"]:
        raise SystemExit("At least one declared validation gate failed; inspect the archived evidence")


if __name__ == "__main__":
    main()
