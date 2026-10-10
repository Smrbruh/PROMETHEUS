# Phase 1.1 reference-label protocol

Defined 2026-10-11 (Asia/Qyzylorda), before new refinement results. This is a separate protocol and evidence namespace. No Phase 1 source, test, report, manifest, failed attempt or output is replaced. The baseline inventory is [preservation_before.json](../results/phase11/preservation_before.json); new files live under new module/script names and `results/phase11/`.

## Reproduction and frozen identity

Run the original six test files (104 tests) and the documented validation, .40/.35/.45 benchmark commands and h1 study into fresh Phase 1.1 directories. Compare scientific outputs, arrays, histories excluding timing, stopping causes, and labels against the originals. Timing and source inventories may differ because new modules exist. Hash every preexisting source/evidence file before and after. Never invoke an old reporting/check script with its default destination if that would replace evidence. New test XML and replay figures have new destinations.

The primary catalogue contains the existing B0/B1/B3_binary/B3_threshold at each volume .35/.40/.45. B0 is the over-budget calibration and B1 the simple comparator. Their historical candidate name, physical-field hash, original source mesh and material/BC definitions are retained. Newly evaluated records additionally identify this phase and their source archive. A reference label is new evidence about that same frozen geometry; it does not rewrite its old status.

Every changed threshold or repaired geometry receives an immutable SHA-256 identity tying source field, geometry, mesh and generation rule together. Exact duplicate geometry may carry an explicit alias, but provenance remains available. No relabeling, refiltering, reranking, smoothing or added support occurs during mesh refinement.

## Numerical refinement and resources

Use the unchanged Q4 and P2 kernels, material, thickness, support and consistent right/down load bases. Load magnitude remains 1,000 N. Binary void cells are removed. First reuse the reproduced h=2/1/.5 histories. For each originally unresolved, connected binary geometry, evaluate h=.25, then .125, then .0625 mm if the previous level still fails and resources permit. Newly repaired geometry starts at h=2,1,.5 and follows the same extension. No extra refinement is required for already accepted B0/B1 unless a new geometry ambiguity is found. Every planned level, success, exception, timeout or resource deferral is recorded; a skipped finest level cannot be treated as a passed solve.

Stop after all declared response settings pass on the last three available consecutive nested levels, or a declared resource/solver limit prevents further progress. Uniform refinement subdivides the exact h2 union of squares; it does not improve the geometry into a smooth manufactured boundary. The study distinguishes linear-system residual convergence from spatial resolution of this fixed mathematical domain.

Resource policy: one numerical worker at a time, numerical backends single-threaded, maximum 600 seconds per combined Q4/P2 mesh attempt, 8 GiB observed child working-set ceiling, and maximum 1.5 million estimated P2 DOFs. Retain the existing 4 GiB conservative P2 COO-assembly preflight; factorization needs additional headroom. Preflight estimates are estimates, not measured memory or condition numbers. A parent monitor can terminate a child at the time/memory ceiling; polling is not a guaranteed OS hard limit. Resource failures remain unresolved. No tolerance is relaxed to finish on this computer.

## Unchanged tolerances and new labels

The original algebraic checks remain: residual <=1e-8, stiffness skew <=1e-12, energy <=1e-8, force/moment balance <=1e-7; Gram reciprocity <=1e-8 before extracting its equivalent symmetric quadratic part. Volume excess is allowed only up to 1e-6 times full domain volume. Successful factorization alone is not acceptance.

For **each** response setting, the final three nested levels must have matching geometry/problem provenance and successful Q4/P2 solves. Last-step relative compliance changes for both methods must be <1%; maximum pad-displacement changes <2%. Finest same-h Q4/P2 differences must satisfy <1% compliance and <2% displacement. Compliance changes must decrease relative to the preceding step, and compliance must not decrease by more than the existing 1e-8 relative roundoff allowance. No monotonicity theorem is assumed for a nodal maximum. Cached differences are recomputed. These are empirical resolution criteria, not rigorous continuum-error bounds; no extrapolated value is substituted for a computed response.

Report separate predicates:

- `numerically_accepted`: the stated fixed linear-model response passed the complete recorded refinement gates.
- `unresolved`: missing levels, failed algebraic/refinement checks, resource limits or an ambiguous domain interpretation prevent the label.
- `geometry_rejected`: disconnected retained material, missing pads or a missing support/load boundary violates the original geometry rule.
- `over_budget`: actual volume violates the budget, regardless of numerical convergence.
- `eligible_reference`: numerical acceptance, valid unambiguous geometry and budget compliance all hold. This is not safety or manufacturing approval.

A globally edge-connected pixel set can still have a vertex with exactly two opposite solid quadrants and two void quadrants. The shared-node mesh then imposes a local point-contact constraint. Count and locate these configurations. If present, retain the original tied-node numerical calculation as a diagnostic but withhold a continuum reference label pending a declared domain/contact interpretation. Never silently split nodes or remove material from the original candidate.

## Reconstruction and additional settings

For every original optimized density field, examine thresholds **.35, .40, .45, .50, .55, .60, .65**, always preserving pads. Archive every resulting mask and identity, actual mass, signed mass error relative to budget and density field, connected components, disconnected area and point contacts. This is a geometry sensitivity sweep; it does not fabricate mechanics results for unevaluated masks.

Feature diagnostics are exact horizontal/vertical maximal contiguous pixel-run distributions and their physical widths, with disconnected areas and local contacts. These reliably describe the raster but are **not** Euclidean medial-axis thickness, minimum load-bearing ligament, a process limit or manufacturability. No smooth fillet radius is inferred from a square grid.

Explore one predefined repair only: in the .35 ranked binary candidate, retain the edge component containing both pads and remove floating components. Do not add material to restore the budget or otherwise change the supported component. Assign a new immutable ID, record removed cells and mass, charge generation and all verification, and make no exact-equal-mass claim if mass changes. If no single component contains both pads, reject this repair rule.

Additional mechanical benchmark settings are load orientations **-20° and +20°**, alongside nominal 0°, each at 1,000 N, on the same frozen candidates. Derive each displacement and compliance from the separately solved right/down basis on every mesh. Apply the same refinement/independent-comparison gates separately to these orientations; their labels do not inherit nominal acceptance. This exact linear reduction consumes no new factorizations; charge postprocessing. These are related loading scenarios, not independent geometry families or evidence of statistical generalization. Keep within-volume rankings separate, show raw ordering separately from eligible ordering, and use the existing 2% reporting-only tie margin without treating it as a certified error bound.

## Large displacement and decision

From archived nodal solutions compute the displacement gradient H, linear strain sym(H), and omitted Green–Lagrange term 0.5 H.T H at Q4 Gauss points. Report area-weighted norms and percentiles globally and outside the existing 6 mm exclusions, plus rotations. This diagnoses kinematic truncation in a **linear solution**; it is not a nonlinear equilibrium solve or a physical failure prediction. Do not infer a safe load or classify physical safety from displacement/length alone. A finite-rotation/geometrically nonlinear plane-stress comparison is the next model if the neglected term is material; actual fixtures/material and out-of-plane effects remain separate questions.

The decision will compare observed resolution/geometry mechanisms and actual costs. Recommend against a learned selective verifier for this catalogue if deterministic validation and refinement suffice and are affordable. If reliable labels cannot be obtained, defer learning rather than training on an arbitrary finest mesh. Neither outcome establishes algorithmic novelty or starts Phase 2.
