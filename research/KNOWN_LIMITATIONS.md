# Phase 1 known limitations

Research date: 2026-10-10. This document separates limitations of the numerical model, the acceptance evidence, and the proposed research claim. The measured values below refer to the canonical 40% run in [its immutable summary](../results/phase1_v040_verified/summary.json), unless stated otherwise. Phase 0 references and problem definitions remain authoritative historical specifications.

## Numerical agreement is not physical validation

The implementation solves a synthetic, static, small-strain, linear-isotropic, plane-stress plate with perfect top-edge fixation. Its material constants are proposed benchmark inputs, not measured or certified properties of a supplied alloy or manufactured part. No physical experiments have been performed. Safety, strength, fatigue life and manufacturing feasibility are not assessed.

The model excludes joint compliance, bolts, contact, plasticity, fracture, fatigue, geometric nonlinearity, impact, material variability and out-of-plane instability. Agreement between meshes and assemblers cannot detect errors shared by these assumptions. The plate is not a validated model of an arbitrary three-dimensional bracket or mounting arrangement.

The nominal finest P2 maximum pad displacements are **13.7412 mm for B1** and **1.94635 mm for B3**. Relative to the 120 mm reference dimension, these are approximately **11.45%** and **1.62%**. Both exceed the protocol's 1% model-review trigger; B1 in particular gives a strong reason to question the small-deformation approximation at the proposed 1,000 N load. This trigger is a provisional review rule, not a theorem defining a safe load or a substitute for a nonlinear model. A numerically accepted linear solution may therefore remain physically inappropriate for the proposed loading.

## Acceptance applies to specified quantities and representations

The canonical complete fixed-field study numerically accepts B0, B1 and relaxed-density B3 under the declared **nominal compliance and maximum pad-displacement** refinement/agreement criteria. B0 remains outside the 40% volume budget. B2 is a coarse grey-density diagnostic and has no complete independent acceptance sequence. Passing `solver_passed` alone establishes only the recorded discrete algebraic checks.

Both **B3_binary and B3_threshold remain `resolution_unverified`**. Their fixed-geometry Q4 compliance changes between h=1 and h=0.5 are approximately 2.90% and 2.98%, above the 1% gate; finest Q4/P2 compliance differences are approximately 1.43% and 1.46%. Their unresolved status must not be replaced with acceptance because they look connected, satisfy a mass constraint, or have a favorable coarse objective. These are measured failures to meet numerical resolution criteria, not demonstrated physical failures.

The full acceptance flag does not independently establish convergence of the angular worst case, unrestricted displacement maximum, stress statistics, or every local field. Maximum displacement on all mesh nodes can be dominated by ersatz regions; the separately reported material-region maximum uses the declared rho>=0.9 region. Compliance, pad displacement, material-region displacement and stress retain different meanings.

Displacement maxima are taken over the stated mesh nodes, not certified continuous-domain maxima. Stress maxima and percentiles are sampled at Q4 integration points. The unrestricted density-model peak includes every cell with positive physical density, including soft ersatz cells; the regular-region statistic uses the stricter rho>=.9 material definition.

The completed separate h=1 optimization-mesh study audits its frozen result at h=1 and h=0.5. It converges in 288 evaluations and changes common finest P2 compliance by -1.5574% versus h2 optimization. Its two levels provide useful comparisons but cannot establish a three-level refinement trend or grant the full h=2/1/0.5 acceptance status. Reoptimization and frozen-field refinement answer different questions. Similar density images are not evidence of mesh-independent optimization. See [the study record](../results/phase1_h1_reoptimization/summary.json).

## Scope of solver independence and correctness checks

The P2 triangle implementation has separate shape functions, constitutive construction, stiffness assembly and load/boundary reconstruction. It still receives the same domain and frozen cell coefficients and shares NumPy, SciPy, sparse direct linear algebra, and some data conventions with the Q4 path. Both implementations were created within this project. The result is useful implementation diversity, not a third-party solver certification or independent material/geometry truth.

Analytical tensile and manufactured-solution fixtures, rigid modes, equilibrium, superposition and derivative checks provide evidence of correctness within their coverage. They are not exhaustive proofs against every coding defect. Q4 energy checks include assembled-matrix energy and the sum of element energies. P2 energy diagnostics use its assembled system; this is not an additional third independent energy implementation or an experimental energy measurement.

LU pivot ratios are recorded as numerical-health diagnostics. **They are not rigorous stiffness condition-number estimates**, and small residuals do not bound forward error in a poorly conditioned model. No general a posteriori error certificate or systematic alternative-ordering sensitivity study is provided. The fixed-field Emin sweep addresses one numerical-model dependence and does not eliminate every near-mechanism or soft-material artifact.

SIMP/filter/OC equations are attributed to established literature and were independently implemented. This phase does not reproduce the complete published top88/MBB benchmark or establish agreement with its full reference optimization history. Agreement on local algebra and the specified L-plate tests should not be described as a published benchmark reproduction.

The finite-difference audit retains well-scaled nonnegative directions and separate sign-changing directions. A first test assertion incorrectly required every sign-changing direction to show decreasing error between the two largest steps, even when cancellation and solve roundoff dominated its tiny truncation error. That assertion was corrected with a mathematical explanation; the actual compliance requirement of relative error below 1e-4 at two adjacent steps was retained. Some smallest-step errors increase. Their complete curves are archived, rather than presenting one selected difference as a proof.

## Geometry, stress and manufactured interpretation

SIMP grey densities are a computational relaxation. They do not describe a qualified isotropic porous material, and void stiffness is a numerical device. Filtering at a 6 mm radius does not guarantee a 6 mm minimum ligament, manufacturability, or tolerance robustness. Binary conversion changes the model and creates a separately identified candidate; mass and performance must be measured again.

Edge connectivity and pad preservation are checked. No certified CAD reconstruction, filleting, corner-radius measurement, minimum-ligament measurement, machining-access analysis, process-tolerance study or manufacturing-cost model has been completed. A union of mesh squares remains a pixel geometry. Connectivity alone does not rule out an impractically thin or fragile connection.

The sharp L corner, clamping endpoints, load-patch endpoints and additional corners in reconstructed shapes can produce stress singularities. Refinement may raise unrestricted peak stress even while compliance stabilizes. The area-weighted 95th-percentile stress outside the specified 6 mm exclusions is a diagnostic only; those exclusions do not cover every new pixel corner and cannot establish strength. No yield allowable, safety factor, stress feasibility constraint or stress-convergence acceptance is claimed.

## Research and computational claims remain limited

The exact angle envelope concerns a fixed, discretized, linear model and the prescribed unchanged load patch over [-20,20] degrees. Its exactness does not validate the mechanics model, establish a probability of failure, or cover moving loads, uncertain supports, contact or nonlinear response. The two-RHS reduction is the required established comparator for this load space; extra sampled angle solves do not create a fair efficiency advantage for an AI selector.

The .35, .40 and .45 volume settings share one geometry, material law, support and load definition. They are related numerical sensitivity cases, not independent samples proving robustness or generalization. Candidate reconstructions and nearby optimization iterates are likewise correlated. The equal-mass selection problem must remain separate within each volume family.

There is no application-derived displacement or strength threshold. Response false-feasibility rates are therefore undefined, not zero. The implemented numerical/geometry predicates can support their own clearly labeled accounting, but cannot stand in for structural safety labels. An unresolved fine reference cannot be used as reliable hidden truth for a learned selector.

No learned selective-verification method, knowledge-gradient implementation, uncertainty-learning comparison or powered stochastic/held-out study has been completed. Two binary alternatives do not establish an allocation-strategy advantage. SIMP is the common candidate generator, not itself a competing allocation policy. PROMETHEUS's algorithmic novelty and superiority remain unproven, and Phase 2 remains gated by the decision memo.

The canonical run measured **94.02 seconds** wall time and approximately **1.99 GiB** peak process working set, including its recorded optimizer repetitions, candidate audits and derivative checks. This is one measured workflow on the recorded machine; it is not a runtime guarantee, a universal hardware requirement, or evidence of a selector's speed advantage. Development time and separately aborted runs are not part of that successful-run wall time. Backend single-thread settings and an 8 GiB checkpoint guard do not impose an operating-system hard memory limit or forecast all sparse-factorization fill.

Bit-identical canonical density arrays and numerical histories were reproduced on this environment. Cross-platform or cross-version bit identity was not established. Runtime can change with machine load, package/backend builds, memory pressure and file I/O.

## Failure evidence is retained

The earlier [aborted run record](../results/phase1_v040/optimization_1/failure.json) preserves a Windows `ctypes` process-handle conversion error in memory measurement. It was an instrumentation failure, not evidence of a successful baseline or a mechanics validation failure. The successful rerun uses corrected handle declarations and its own source snapshot. Do not combine its successful-run timing with an invented duration for the earlier abort.

The canonical successful-run summary records zero failed audit attempts, while the two reconstruction resolution failures remain explicit verification outcomes. The experiment runner retains known completed primary results and arrays if a later independent stage fails, records failed attempts and their available cost, and does not invent completed factorization counts for unknown stages. Failure ledgers and unresolved statuses are part of the evidence and must accompany any selected results.
