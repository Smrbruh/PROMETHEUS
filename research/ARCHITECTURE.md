# Recommended research architecture
Date: 2026-10-10. This is a design proposal, not an implemented application.

## Decision and environment
Use **Python orchestration and numerical arrays with existing compiled solver libraries**, starting on CPU. Implement only the benchmark and its verification harness after Phase 0. Neither a Rust rewrite nor a custom C++ FEM kernel has a demonstrated need.

Inspection found an empty directory, no Git repository, no source, project license, dependency manifests, tests or prior implemented functionality. No applicable ancestor AGENTS.md was found. Windows 11, Intel N100 (4 cores/4 threads), 15.75 GiB RAM and approximately 81.32 GiB free disk were observed. Python 3.14.4 is installed with NumPy 2.4.6, SciPy 1.17.1, pytest 9.1.0 and matplotlib 3.10.9. Torch is discoverable but was not exercised; Gmsh, scikit-fem, FEniCS, pymoo and BoTorch were not discoverable in that interpreter. Executables for Git, CMake, Rust and Cargo exist; compiler toolchains and GPU suitability were not validated. These facts do not establish binary compatibility of future dependencies.

Create an isolated environment when implementation starts. Select the newest mutually supported Python release after testing required wheels; do not assume every library supports the installed 3.14 runtime. Pin resolved versions, platform, BLAS backend and solver options. Do not modify the global installation to make the architecture appear ready.

## Language and technology choices
| Concern | First choice and rationale | Alternative / gate |
|---|---|---|
| Numerical arrays/sparse algebra | NumPy/SciPy; installed, compiled kernels, inspectable residuals | Profile before custom kernels; direct sparse factorization first |
| Reference optimization | Independently implemented SIMP + density filter + OC equations; compare DTU algorithm | DTU educational source is a reference, not presumed redistributable code |
| Verification FEM | Separately assembled second-order triangular plane-stress solution using scikit-fem | CalculiX separate executable for later 3D; FEniCSx if PDE generality justifies Windows deployment effort |
| Meshing | Structured Q4 for optimization; conforming triangular mesh for verification | Gmsh Python API/CLI once compatibility and license reviewed |
| Geometry | Arrays + explicit masks initially; immutable contours for reconstructed geometry | Open CASCADE through FreeCAD/CadQuery only when solid CAD is required |
| Gradient/constrained optimization | Analytic adjoint derivative and OC for v0 | SciPy constrained methods for small parameter problems; audited MMA implementation for many constraints |
| Pareto/discrete design | Fixed-budget volume sweep first | pymoo for small-dimensional mixed/discrete study, count all evaluations |
| Uncertainty | Exact load-space reduction, deterministic scenario sweeps | Small GP/active reliability only after discrepancy pilot; BoTorch optional |
| Visualization | matplotlib diagnostic figures, optional mesh viewer | Display is a consumer of saved results, never the solver |
| Tracking | Human-readable manifest + NPZ/CSV/JSON and immutable run directories | Database/MLflow/distribution only after measured coordination need |
| Testing | pytest numerical invariants and independent reference outputs | No UI or extensive infrastructure tests before any UI exists |
| CAD interchange | DXF for planar contours; STEP for real solid model later | STL is tessellation, not exact CAD; verify watertightness and units separately |
| Cost/manufacturability | Process predicates and transparent parameterized cost breakdown | Supplier/process calibration before currency estimates |

See [COMPETITOR_ANALYSIS.md](COMPETITOR_ANALYSIS.md) for official capability and license evidence. Expected software choices: permissively licensed NumPy/SciPy/scikit-fem; GPL implications for Gmsh/CalculiX; LGPL-family obligations for FreeCAD/Open CASCADE; confirm exact versions and distribution mode before bundling. Calling a separately installed tool does not erase obligations or establish a legal conclusion. Commercial tools are optional documented comparators, not prerequisites.

Rust offers memory safety and good typed interfaces but adds bindings and a smaller ready-made research ecosystem. C++ offers mature geometry/solver kernels already available through bindings, but a new solver creates verification burden. Use either only after profiling identifies a concrete bottleneck or an existing reusable library requires it.

## Module boundaries and contracts
All records carry schema_version, units, immutable content hash, creation provenance and explicit status. Plain Python functions/dataclasses are sufficient. No service API is implied.

| Module | Contract | Required checks |
|---|---|---|
| Problem specification | ProblemSpec → domain, material, load cases, constraints, prescribed displacements | Positive E/thickness; valid ν; units; supported BCs; required fields |
| Candidate generation | ProblemSpec + GeneratorConfig + seed → Candidate | Representation, parameters and parent history; no performance claims |
| Material model | MaterialSpec + physical density → constitutive tensor | State plane stress/strain/3D; no undocumented alloy allowables |
| Geometry transform | Candidate + TransformConfig → GeometryArtifact | Mapping density→contour→solid; areas/masses before/after, islands, holes |
| Meshing | GeometryArtifact + MeshSpec → MeshArtifact | Element orientation/Jacobians, boundary tags, scale, quality and hash |
| Physics solve | MeshArtifact + ProblemSpec + SolverSpec → SimulationResult | Success/error status; residuals, reactions, energy, runtime and fields |
| Optimization | Candidate + objective/sensitivity interface → OptimizationHistory | Budget charged, physical volume, gradients, stopping cause |
| Cheap constraints | Candidate/GeometryArtifact → PredicateReport | Connectivity, pads, thickness/process predicates; unsupported checks explicit |
| Independent verification | Frozen Candidate + VerificationSpec → VerificationReport | Separate assembler/element/mesh; unchanged loads and compared geometry semantics |
| Uncertainty analysis | ProblemSpec + uncertainty set + solved data → ResponseEnvelope | Exact reductions where applicable; epistemic vs probabilistic assumptions |
| Adaptive selection | Query history + budget → QueryProposal | Proposal only; cannot change saved ground truth or relax constraints |
| Benchmark runner | Frozen manifest + algorithm spec → RunRecord | Seeds, held-out split, failure policy, resource accounting |
| Analysis/export | Saved RunRecords → tables/plots/CAD | Reads authoritative records, includes invalid/abstain outcomes |

A solver failure returns status such as invalid_geometry, singular_system, nonconverged, resource_exhausted or unsupported_model. It must not return zero compliance, silently clip stresses or masquerade as a feasible result. Optimizers may use a declared rejection/penalty policy, while raw failure remains visible.

## Data flow
ProblemSpec → generation/optimization → frozen candidate → geometry/mass/cheap checks → independent mesh and solve → verification report → comparison.
Optional adaptive allocation chooses which frozen candidate to send through verification next. Visualization and language-model summaries read the reports. They cannot write approved numerical fields.

Optimization analysis and independent verification may share units/schema but must not share the same stiffness assembly kernel or claim independence from merely changing a plotting resolution. First compare identical density fields; later compare an explicit reconstruction. Keep numerical, geometric and physical discrepancies separately labeled.

## Minimum candidate provenance
A run directory should eventually contain:
- manifest.json: problem ID/revision, hashes, timestamps, command, machine/OS, dependency lock and thread settings.
- specification.json: SI-compatible unit convention, geometry, material source/status, load and BC coordinates, uncertainty model, objective and constraints.
- design.npz plus geometry: raw design variables, filtered physical density, pads, reconstruction parameters, parent/seed.
- mesh and boundary tags with mesher version/options/hash.
- history.csv: objective, volume, gradient diagnostics, solver metrics, iteration counts, timings and stopping reason.
- analysis.json and verification.json: fields/metrics locations, tolerance outcomes, explicit feasibility/abstention, warnings and failure records.
- licence/source ledger: algorithm references and third-party version/license notices.

Keep floating-point arrays in binary format with documented axis ordering; never use screenshots as results. Hash canonical files, record non-deterministic software behavior, preserve stdout/stderr and a rerun command. Pinning software helps repeatability but does not replace independent verification.

## PC resource policy
v0 optimization has 1,600 Q4 elements. Primary Q4 refinement at 6,400 and 25,600 cells is staged. The independent P2 triangular mesh is larger: splitting the finest full-domain grid gives 51,200 triangles and roughly 206,000 displacement degrees of freedom before constraints, not the Q4 count. Actual peak memory depends on ordering and factorization fill, not just DOF count. Measure one mesh at a time before the full study. A provisional process-memory cap of 8 GiB leaves operating-system headroom; stop cleanly before paging dominates and log the cap. It is a conservative operating policy, not a measured performance bound.

Use one benchmark process at a time, fix BLAS/OpenMP threads, save cold/warm timings separately and exclude neither meshing nor surrogate fitting. If fine verification exceeds resources, refine locally or reduce the experimental scope; do not quietly lower fidelity for only one method.

## AI strategy and adoption gates
| Application | Conventional comparator | Evidence required before adoption |
|---|---|---|
| Natural-language requirements | Typed form/schema and engineer review | Constraint extraction error and silent omission rate; valid schema does not establish correct intent |
| Candidate suggestions | SIMP/OC, parameter sweep, simple geometric references | More verified feasible value at equal total budget |
| Verification ordering | Coarse-rank, exact physics, random/space-filling, U-learning, knowledge gradient | H2 held-out selection quality/cost and calibration |
| Surrogate mechanics | FEA, ROM/POD, simple response surface | Error at constraints and OOD cases; data-generation/training cost included |
| Solver anomaly detection | Residual/energy/reaction/mesh checks | Added true detections with bounded false alarms, never substitute for checks |
| Counterexample proposal | Exact load-space oracle and deterministic sweeps | Better found violations per total cost on nontrivial uncertainty only |
| Report drafting | Saved tables and scripts | Every numeric sentence traceable to result hash; no fabricated runs |

AI outputs are proposals. Parsed constraints require dimensional, range and semantic validation before becoming a ProblemSpec. A language model may not invent material allowables, authorize manufacturing or testing, alter acceptance thresholds, declare a failed solve successful, or produce authoritative stresses from text. No neural training is recommended in the first milestone. A GP is optional for H2, with standard non-AI policies evaluated alongside it.

The first architecture deliverable after Phase 0 is a small command-line experiment and its records. The product's eventual UI, services, database, language-model provider and GPU stack are intentionally undecided because the research has not demonstrated a need.
