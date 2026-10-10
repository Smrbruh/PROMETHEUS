# PROMETHEUS-01 baseline experiment

Status: **preregisterable proposed protocol; no bracket optimization result has been produced in Phase 0**. Canonical identifier: `PROMETHEUS-01-v0`. All geometry, units, constitutive equations, filter derivatives, and accounting follow [MATHEMATICAL_FOUNDATION.md](MATHEMATICAL_FOUNDATION.md). These are deliberately proposed settings, not a claim that this is an established public benchmark.

## 1. What the first experiment establishes

The first experiment should establish whether a small, independently checked SIMP pipeline can generate and compare meaningful plate candidates on the available PC, and whether a geometric conversion creates measurable verification discrepancies worth researching. It cannot establish general bracket safety, novelty, or an advantage of AI.

Record the full solid solution first; then an equal-volume simple structure; then a deterministic density-filtered SIMP/OC solution. The algorithmic reference is the authors' educational top88 density-filter approach [M01, [official reference](https://www.topopt.mek.dtu.dk/apps-and-software/efficient-topology-optimization-in-matlab)]. Its source footer reserves rights [T46, [author-hosted source](https://www.topopt.mek.dtu.dk/-/media/subsites/topopt/apps/dokumenter-og-filer-til-apps/top88.m)]; do not vendor it as permissively licensed code. Independently implement the published equations with attribution. Exact pixel-for-pixel agreement with its MBB example is a useful regression only when all of that example's settings match; it is not independent physical validation of this L plate.

## 2. Frozen input manifest

| Field | v0 value |
|---|---|
| Geometry | L plate: 120×80 mm bounding box, 40 mm legs; 6,400 mm² active area |
| Thickness | 5 mm |
| Preserved pads | [0,40]×[70,80] and [110,120]×[10,30] mm² |
| Restraint | ux=uy=0 only on y=80, x∈[0,40] |
| Load | Uniform on x=120, y∈[10,30]; total 1,000 N downward |
| Material | Synthetic E=70,000 MPa, ν=0.30, density 2.70×10⁻⁶ kg/mm³ |
| Volume | 0.40 of total active domain, including 600 mm² of preserved solid pads |
| Element | Full-integration plane-stress Q4, h=2 mm, 1,600 active cells |
| Optimizer | OC; move=0.20; x∈[10⁻³,1]; p=3 throughout |
| Filter | Density filter, r=6 mm; pad-aware normalization defined in mathematical document |
| Projection | None, β=0 identity; no continuation |
| Void stiffness | Emin/E0=10⁻⁶ |
| Solver | Double precision, sparse direct solve, explicit free-DOF elimination; pin actual package/build/version in execution manifest |
| Baseline budget | At most 300 evaluated designs per optimization run, including the initial design |
| Initial field | Constant free x=α, pad x=1; bisection on α to match filtered physical volume |
| Default randomness | None; deterministic geometry ordering and matrix assembly |
| Angular audit | Exact θ∈[−20°,20°] compliance envelope from two RHS; nominal θ=0 |

Use SciPy's sparse LU path initially; record ordering, actual backend, thread count and numerical library configuration. Do not silently switch solver backends across methods. A factorization can be reused for multiple RHS for a frozen design. Sparse direct cost and memory depend on fill-in, so no constant seconds-per-design claim is made before timing. The h=2/h=1/h=0.5 domains contain 1,600/6,400/25,600 Q4 cells, respectively.

The base domain uses rows from low y to high y and columns from low x to high x; only retain cells whose centres lie in Ω. All specified pad and traction endpoints align with these three grids. A cell belongs to a pad when its centre lies strictly inside that pad or on its boundary according to an explicitly recorded tolerance; no centre lands on a pad boundary for these grids. Node IDs need not match another solver, but exported coordinates, connectivity and boundary tags must.

## 3. Reference designs and fair comparison

| ID | Construction | What its comparison means |
|---|---|---|
| B0 | Full solid Ω, all ρ=1 | Debugging and stiffness reference; violates 40% volume budget and cannot win the equal-mass competition. |
| B1 | Binary thick L-path, including preserved pads, at exactly 40% cell volume | Simple geometry comparator; construction below avoids arbitrary manual tuning. |
| B2 | Constant raw design α initialized to satisfy filtered volume, preserved pads solid | Numerical grey-density diagnostic, not a fabricated part or optimized comparator. |
| B3 | Density-filtered SIMP/OC from B2 | Established optimization baseline; report relaxed and converted binary outcomes separately. |

For B1, define the polyline through (20,75), (20,20), and (115,20) mm. Start with all 150 pad cells solid. Sort the 1,450 designable h=2 cells by Euclidean centroid distance to this polyline, breaking ties by increasing y-index then increasing x-index. Fill the first 490 cells. Thus 640/1,600 cells are solid and the total volume is exactly 12,800 mm³. Freeze this h=2 union of squares when refining its analysis mesh. Do not regenerate a different geometry during solution convergence. This deliberately simple baseline may have inferior load transfer; it is a comparator, not an engineer-approved bracket.

For equal-mass **binary** reporting, convert a B3 physical-density field into solid cells by retaining both pads and the 490 highest-density designable cells, using the same lexicographic tie-break. Archive the exact selected indices. This is a declared conversion rule, not optimization or smoothing. Also report the conventional ρ≥0.5 threshold at its **actual** volume as a sensitivity diagnostic. No silent geometry repair or volume rescaling is allowed. A disconnected or point-contact-only result is a failed geometric candidate, even if the ersatz solve returns finite compliance.

For density comparisons, evaluate B2/B3 with the same SIMP interpolation and actual physical density. For binary comparisons, evaluate B1 and the converted B3 with material only in solid cells and traction/support geometry preserved. Do not claim that lower relaxed compliance automatically transfers to a manufactured part. If a repaired or smoothed shape is considered later, it receives a new candidate ID, measured area, remesh and complete re-verification.

## 4. Iteration, stopping, and failure

At each evaluated design: filter; override pads; check volume; assemble and solve; compute objective and sensitivities; archive metrics; then decide whether to stop or perform an OC update. Return the last **evaluated** candidate, never an unanalysed update.

After at least 30 evaluated designs, the proposed convergence flag requires both conditions to hold over the last ten transitions:

1. Maximum raw-density change ≤0.01 at each transition.
2. Relative spread of the last eleven compliance values, (max C−min C)/max(|Ccurrent|, 10⁻¹² N·mm), ≤10⁻³.

Volume violation must be ≤10⁻⁶ VΩ and all solve checks must pass. The bisection itself targets volume error ≤10⁻⁸ VΩ or relative multiplier-bracket width ≤10⁻⁸, with at most 100 bisections, accepting the feasible endpoint. A stop at 300 evaluations without convergence is `budget_exhausted`, not `converged`; report objective, remaining change and best feasible evaluated incumbent. If a feasible move-limit candidate cannot be found, mark the iteration failed rather than changing the budget.

These optimization stopping choices are **provisional numerical protocol settings**. Pilot calibration may change them once, before comparisons are frozen; the revision requires a new manifest and rerunning every affected comparator. Strong residual, conservation, and derivative checks in the validation protocol are implementation requirements, not tuning knobs for improved results. Terminate and retain diagnostic evidence on nonfinite data, failed factorization, invalid Jacobian, unsupported load patch, or residual/conservation failure.

## 5. Outputs and interpretation

Save a row per evaluated design and a final report containing:

- Nominal compliance in N·mm and C/Csolid on the same analysis definition; relaxed and binary values in separate columns.
- Total physical volume, total mass, preserved volume, variable volume and volume violation.
- Load-weighted displacement C/F in the nominal load direction, average load-patch displacement, maximum pad displacement, and maximum displacement on the declared material region.
- Exact angular worst-case compliance and its angle; three fixed report cases θ=−20°,0°,20° are illustrative summaries, not the complete angular maximization.
- Unrestricted peak and explicitly labelled regular-region stress diagnostic, plus their mesh trends; no stress safety factor.
- Grey indicator G=ΣD 4ρ(1−ρ)Ae/ΣD Ae; binary topology/connectivity status; conversion-induced area and compliance changes.
- Factorizations, RHS solves, derivative evaluations, geometry operations, iterations, wall time, CPU/thread settings, peak memory and all failures.
- Mesh/solver/transfer checks and their tolerances; no single composite score that hides failed constraints.

The load-weighted displacement C/F is a scalar associated with the prescribed distributed loading, not necessarily the maximum nodal displacement. Count all preprocessing, candidate generation, surrogate fitting and failed evaluations when measuring end-to-end algorithm cost. Plot objective versus both elapsed time and factorization/RHS counts; reusing a factorization is a real advantage that the accounting must preserve.

## 6. Independent verification and mesh dependence

First freeze B0, B1 and the final B3 density field at h=2. Subdivide each parent cell into h=1 and h=0.5 children while retaining the parent's physical density; **do not refilter, threshold again, or reoptimize**. This checks the solution for a fixed material coefficient field. The independent verifier assembles plane-stress elasticity with quadratic triangles, respects all cell-material interfaces, integrates the same line traction and thickness, and uses a separate assembly code path. scikit-fem is a plausible lightweight verifier, with an open software publication and elasticity examples [M05, [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.02369); M06, [official gallery](https://scikit-fem.readthedocs.io/en/latest/listofexamples.html)]. The implementation and its plane-stress convention still require verification.

For binary geometry, remove void cells entirely in the verifier. Split each retained parent square consistently into triangles, refine conformingly, and retain physical boundary tags. If disconnected free components cause rigid modes, fail the geometry; do not ground them or add weak springs. Higher-order triangles on the same shape provide different approximation and assembly, but not independent material or boundary-condition truth.

Separately reoptimize on h=1, preserving r=6 mm, volume fraction and all geometric dimensions, to test **optimization mesh dependence**. If needed, h=0.5 reoptimization is a later experiment after memory/time measurement. Compare all optimized outcomes with common verified metrics, not only their own mesh objectives. A geometrically similar image or a fixed radius does not prove mesh independence. [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md) defines numerical gates.

## 7. Repetition, holdout and research-budget discipline

The deterministic baseline first requires repeatability, not artificial statistical replication. Run the identical manifest twice to verify reproduced values and hashes. Later method comparisons use the same candidate catalog, train/test split, initial observations, solver versions, geometry conversion, and budgets for all methods. They may not discard failures or let the proposed method inspect held-out verification labels.

If random initialization or stochastic selection is introduced, use the common development seeds `[0,1,2,3,4,5,6,7,8,9]`. Store the actual initial arrays and generator specification rather than relying only on seed portability. These ten runs are a feasibility pilot, not a powered significance claim. Choose and freeze the confirmatory repetition count from the pilot's variance and the engineering-relevant effect, before unblinding test outcomes. A new held-out geometry family, attachment model or load-patch location tests transfer; more angles on the same frozen linear plate do not provide independent expensive physics tasks.

The 300-design budget applies to the elementary B3 optimizer. High-fidelity verification and adaptive-selection budgets for the proposed research hypothesis must be calibrated from the pilot, separately preregistered, and inclusive of warm starts/training. No numerical percentage improvement or time saving is asserted in advance. See [RESEARCH_HYPOTHESES.md](RESEARCH_HYPOTHESES.md) for the selection study and [IMPLEMENTATION_ROADMAP.md](IMPLEMENTATION_ROADMAP.md) for gates before building it.

## 8. Conditions that change the recommendation

Stop and repair the numerical pipeline if physics and derivative checks fail. Narrow the hypothesis if binary conversion dominates every difference or if all discrepancies are captured by a simple deterministic indicator. Do not train an expensive selector if independent solves are cheap enough to exhaust the finite catalog. Proceed to an adaptive-verification experiment only if a frozen discrepancy audit establishes nontrivial, reproducible uncertainty in the relevant engineering decisions and a useful cost trade-off against simple baselines.
