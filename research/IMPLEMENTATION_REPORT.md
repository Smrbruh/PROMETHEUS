# Phase 1 implementation report

Research date: 2026-10-10. PROMETHEUS-01 implements the Phase 0 mathematical definition as a CPU Python numerical benchmark. Established SIMP, density filtering and optimality criteria are independently implemented from the documented equations with attribution to Andreassen et al. (2011), DOI 10.1007/s00158-010-0594-7, already recorded as M01. No top88 source was copied. Novelty remains unproven.

## Audit and implementation sequence

The audit read the decision/status/charter, hypotheses, mathematical foundation, baseline, architecture, validation protocol, roadmap and references before selecting the implementation. The package supplied an unusually specific problem definition: no replacement geometry or easier load case was needed. The plan was to (1) freeze units and interfaces, (2) verify a Q4 elasticity kernel and the existing envelope mathematics, (3) implement/filter/check SIMP sensitivities and OC, (4) build a separately assembled P2 reference, (5) run deterministic candidates and fixed-field refinement, and (6) report evidence before deciding on selective verification. No Phase 2 method was implemented.

## Frozen problem and explicit choices

The domain is ([0,120] x [0,40]) union ([0,40] x [40,80]) mm, area 6,400 mm². Thickness is 5 mm. Synthetic E=70,000 MPa, nu=.30 and density=2.70e-6 kg/mm³ imply full mass 86.4 g. Both displacement components are fixed only on the top edge y=80, 0<=x<=40. Uniform line traction on x=120, 10<=y<=30 integrates to 1,000 N downward, with origin moment -120,000 N mm. Pads [0,40] x [70,80] and [110,120] x [10,30] remain solid. The 40% budget includes their 3,000 mm³: total 12,800 mm³, mass 34.56 g.

The main Q4 mesh has h=2 mm, 1,600 cells. Plane stress uses full 2x2 Gauss integration and engineering shear strain. Sparse CSC stiffness is solved by SciPy SuperLU with COLAMD ordering and explicit essential-DOF elimination. No stiffness averaging, weak springs or silent grounding is used. The regular grid and binary subsets preserve actual boundary segments; a missing load patch is an error.

Density filtering uses active-domain centroids, physical radius 6 mm, area weights and the exact normalization/pad override in the specification. SIMP p=3, Emin/E=1e-6, raw lower bound .001 and OC move .20 remain fixed; no projection or continuation. Constant free raw initialization is determined by physical-volume bisection, not by setting raw density equal to the target fraction. At 40% its value is approximately .329920670483. The optimizer returns the last evaluated design and separately archives its best evaluated incumbent. The 300-evaluation cap includes initialization. Convergence requires at least 30 designs, ten raw-change transitions <=.01, last-eleven compliance spread <=.001, and feasible volume. Exhaustion is reported as exhaustion.

Binary geometry uses edge connectivity and requires every retained component and both pads to connect to the supported component. Point contact is invalid. B1 follows the frozen distance-to-polyline rule. B3_binary retains pads and highest-ranked physical-density cells at exact budget; ties are ordered by y then x. B3_threshold retains rho>=.5 at its actual mass. Void cells are removed for binary solves. There is no repair or smoothing.

## Resolved ambiguities and implementation boundaries

| Issue | Explicit resolution |
|---|---|
| Independent solver package | A separately written quadratic-triangle implementation uses the existing NumPy/SciPy stack. The Phase 0 protocol allowed another audited independent solver; no scikit-fem installation is needed. |
| Angle convention | Basis columns are RIGHT and DOWN with full force magnitude included. Coefficients are [sin(theta),cos(theta)]. This is algebraically equivalent to the Phase 0 upward basis with a negative vertical coefficient; interval remains [-20,20] degrees. |
| Numerical reciprocity versus exact envelope | The algebra utility retains its 1e-12 symmetry gate. The mechanics adapter records raw Gram skew, rejects reciprocity error >1e-8, then passes (Q+Q.T)/2. For any real vector a, a.T Q a equals a.T (Q+Q.T)/2 a exactly. This handles measured solve roundoff without averaging stiffness or changing the physical model. |
| Refinement identity | Children inherit the original physical coefficient exactly. No refiltering, reranking or reoptimization occurs inside the h2/1/.5 audit. Different representations and optimizer grids are separate candidates. |
| Numerical acceptance | Complete provenance, preserved pads, full mesh sequence, algebraic checks, decreasing refinement changes, nondecreasing compliance within 1e-8 roundoff allowance, and finest independent agreement are required. Cached relative gaps are recalculated from raw metrics. |
| Feasibility | Volume/pads/geometry/algebraic checks define optimization feasibility. No displacement or strength limit exists, so response feasibility is undefined. Independent numerical acceptance does not require B0 to meet a mass budget; B0 remains excluded from equal-mass selection. |
| Stress | Gauss-point von Mises values are diagnostics. The regular region excludes physical 6 mm neighborhoods of the five specified singular points and requires rho>=.9. No allowable, safety factor or failure prediction is provided. |
| Resource policy | Single-thread numerical backends, sequential benchmarks, measured 8 GiB process-peak checkpoints and a conservative P2 assembly guard. These are not OS-enforced hard memory limits. |
| Additional settings | Volume .35/.45 were declared before observing results. They are sensitivity cases of the same family, not independent validation populations. |

## Interfaces

| Module | Contract |
|---|---|
| `models.py` | Validated unit-bearing Material, Benchmark and LoadCase; explicit NumericalError |
| `geometry.py` | Mesh coordinates/connectivity/indexing; domain and subsets; physical pads, fixed DOFs and unit load basis |
| `fem.py` | Q4 element/assembly, multiple-RHS solve, reactions/residual/energy diagnostics and stress outputs |
| `optimization.py` | DensityFilter physical field and transpose chain rule; volume gradient; SIMP sensitivity/OC; history and binary geometry predicates |
| `envelope.py` | Exact endpoint/stationary-point maximum of a frozen 2x2 symmetric compliance quadratic over a bounded angle interval |
| `verification.py` | Conforming P2 triangles, separate constitutive/shape/stiffness/traction assembly, constraints, solve and manufactured-solution reference |
| `evaluation.py` | Immutable Candidate; frozen coefficient transfer; common metrics/status/provenance; full-audit acceptance |
| `run.py` | Manifest, snapshots, complete run/repetition/audit workflow, evidence preservation, cost and failure records |

The verifier shares the supplied geometry, coefficient field and sparse algebra backend with Q4; it does not import Q4 element, assembly, constitutive or traction code. This is implementation diversity, not an external certified reference. Smooth analytical fixtures exercise each implementation separately.

## Reproducibility and evidence

Each successful experiment archives source and specification snapshots, SHA-256 hashes, package/backend versions, thread settings, initial/raw/physical/best density arrays, mesh and boundary data, full iteration histories, candidate coordinates/connectivity/forces/displacements, Gauss stress arrays, independent fields, stopping cause, numerical diagnostics and timing. Design hash and problem/source-geometry identity prevent unrelated audit records being combined. Existing output directories cannot be overwritten.

The canonical run repeats optimization from the same archived initialization and compares raw/physical arrays and numerical history after excluding timing. Deterministic identity is a measurement on this machine, not a portability guarantee. Recorded wall time includes evidence writes within the runner; solver timing and factorization/RHS counts are reported separately. A failed independent call retains any completed primary result and its known cost; unsuccessful stages are not counted as invented successful factorizations.

Reproduction commands are in the root [README](../README.md). Actual measurements, failures and remaining limits are in [BASELINE_RESULTS.md](BASELINE_RESULTS.md), [NUMERICAL_VALIDATION.md](NUMERICAL_VALIDATION.md) and [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md). The original Phase 0 documents and references remain historical evidence.
