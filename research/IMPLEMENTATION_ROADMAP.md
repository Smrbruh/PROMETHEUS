# Implementation roadmap and decision gates
Date: 2026-10-10. All milestones below are proposed; Phase 0 did not implement the bracket solver.

## Sequence
| Milestone | Smallest output | Evidence to exit | Stop/revise condition |
|---|---|---|---|
| 0 — This evidence package | Charter, reviewed prior art, mathematical benchmark, tool matrix and adversarial review | Claims mapped to sources; missing evidence explicit | Broad invention claim remains unsupported |
| 1 — Mechanics kernel and checks | CPU plane-stress assembly and unit/energy/patch/gradient checks | [Validation tolerances](VALIDATION_PROTOCOL.md), no unresolved solver errors | Fix mechanics before optimization |
| 2 — Reproducible established baseline | PROMETHEUS-01-v0 SIMP/OC + solid and simple fixed-mass references | Convergence status, exact mass accounting, finite loads, archived settings | Different constraints invalidate comparison |
| 3 — Independent discrepancy audit | Same-field independent solution, then frozen reconstruction with mass/geometry checks | Separately quantified mesh and geometry effects; measured PC costs | No meaningful discrepancy → ordinary verification suffices |
| 4 — H2 preregistered pilot | Small catalogue, exhaustive hidden audit data, simple allocation comparators | Calibrated loss/thresholds, resource budget, sample-size justification | Surrogate overhead or weak correlation defeats premise |
| 5 — Confirmatory method comparison | Frozen closest-method implementations and held-out families | Paired uncertainty intervals, full cost accounting, failures/abstentions | No useful difference against strong baselines → report null and stop method branch |
| 6 — Process-realistic 3D model | Solid CAD, fillets, real supports/joints, material/process parameters | Independent 3D numerical verification and process review | Plane-stress assumptions fail → revise problem |
| 7 — Supervised physical study | Calibrated fixture/sensors/coupons and blinded specimen predictions | Independent measurements and test-specific approval by responsible lab | No safe apparatus, unidentifiable model or mismatch → no autonomous testing |
| 8 — Research prototype scope | Only modules shown useful by evidence | Reproduced benefit and a specific user need | Do not build general product without this evidence |

No calendar duration is promised: Windows dependency availability, solver performance and laboratory access have not been tested. Milestones are acceptance gates, not a commitment to build everything.

## The single next action
Implement and run **one established SIMP bracket baseline with a paired independent verification audit**, first comparing the same density field and then a declared reconstruction. Establish how much of any discrepancy comes from mesh, geometry or model. This has greater information value than training a model: it can eliminate the proposed research need.

Before optimizing the bracket, perform the milestone-1 numerical checks; they are part of the baseline experiment, not a separate product feature. Record a timing/memory pilot and reduce scope if the N100 machine cannot support the specified fine mesh.

## Research implementation order
1. Freeze canonical benchmark specification and input serialization.
2. Verify elasticity and sensitivity calculations against analytic/independent checks.
3. Implement filtering/OC from equations with source attribution and license review; no copied unlicensed reference code.
4. Reproduce a published educational benchmark to check algorithm behavior, then the new L-plate configuration. Matching a different domain is not expected to reproduce a paper's objective number.
5. Add independent geometry/mesh/solver verification and exact angle-envelope analysis.
6. Audit a fixed candidate catalogue; decide whether H2 has sufficient effect and cost headroom.
7. Only then evaluate surrogate-driven allocation and the strongest applicable published comparators.

## Confirmatory protocol before method development
Freeze generator settings, candidate family partitions, all physical/numerical constraints, fidelity definitions, cost limits, hyperparameter tuning allowance, evaluation metrics, seeds and stopping rules in a manifest. Reserve a separate pilot set for thresholds and effect calibration. Record every algorithm attempted, including failed or discarded variants.

Use matched wall-time/CPU-budget plots and counts of factorizations/right-hand sides. Independent whole-protocol verification consumes resources too; show development, training, evaluation and final-audit cost separately. No “speedup” based solely on skipping final validation or hiding training time.

If future industrial cost requirements appear, incorporate them explicitly and validate the cost estimator. Until then “lightweight” and “cheaper” must not be equated.

## Outcomes that justify continued work
A baseline that passes verification; repeatable, meaningful reconstruction/model discrepancy; a measurable budget advantage over strong conventional methods on unseen families; and a precise prior-art delta or useful reproducibility contribution.

If discrepancies are negligible, adopt established optimization and ordinary verification as the engineering solution. If discrepancies are large but unlearnable, increase deterministic verification or improve the representation. If off-the-shelf methods already solve allocation adequately, integrate them and abandon algorithm novelty. If practical loading/manufacturing requirements remain undefined, keep the output a benchmark rather than a product claim.

