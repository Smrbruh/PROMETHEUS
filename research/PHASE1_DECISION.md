# Phase 1 decision: retain the baseline; do not start Phase 2

Decision date: 2026-10-11, Asia/Qyzylorda. Study definition began on 2026-10-10; archived run timestamps are UTC. **NARROW the next research step to numerical resolution and representation. The evidence does not yet justify implementing the selective-verification hypothesis.** Novelty remains unproven.

The project now has a small, reproducible computational baseline: Q4 plane-stress elasticity, an established density-filtered SIMP/OC optimizer, exact two-load angular compliance, a common candidate evaluator, and a separately assembled quadratic-triangle verifier. Analytical, manufactured-solution, sensitivity, failure and provenance checks support these implementations within their stated scope. This is useful numerical research infrastructure, not a validated engineering product.

## Evidence that changes the decision

At 40% volume, optimization converged in 114 evaluated designs, and its arrays and numerical history reproduced exactly. The frozen density field passes the complete nominal compliance/pad-displacement audit. Its Q4 h2 compliance is 1,802.453 N mm; finest P2 is 1,872.108 N mm. The full solid and simple binary L-path also pass the declared resolution gates.

The exact-mass binary reconstruction looks much better in the coarse model: 1,123.267 N mm versus the simple L-path's 12,658.591. Its finest P2 value rises to 1,248.481, but it **still fails** the resolution rules: 2.9017% final Q4 compliance change, 1.1573% P2 change and 1.4274% finest inter-solver difference. The threshold reconstruction also fails. Their raw ordering does not reverse; their eligibility changes. An audited-only numerical selection retains B1 rather than promoting an unresolved response to a trusted result.

The predefined nearby settings expose additional limits. At 35% both reconstructions contain disconnected material and are rejected before solving; the simple L-path misses its refinement gate, so no equal-mass binary design is accepted. At 45% the ranked reconstruction again misses resolution gates, and the threshold reconstruction also exceeds the volume budget (.4525 versus .45). Every density B3 field passes its own nominal-response audit, but that does not transfer to its reconstructed binary geometry.

This supports **ordinary independent verification and explicit representation tracking**. It does not yet support a learned allocation method. There are only two equal-mass binary alternatives per setting; the settings share one model family. Their reference labels are partly unresolved, geometry failures are detected deterministically, and both canonical binary audits together took about 8.49 s in evaluation-stage wall time. No knowledge-gradient, expected-improvement or proposed selective method has been implemented or compared, so no advantage can be inferred.

Numerical acceptance also does not validate the loading model. B1's finest canonical pad displacement is 13.74 mm, about 11.45% of the 120 mm scale, making its small-deformation interpretation especially questionable. No candidate is recommended here for manufacture or real loading.

## Completion evidence

The additional h=1 optimization converged in 288 evaluations and reduced common finest P2 compliance by 1.5574% versus h2 optimization. Its two-level diagnostic audit passes its thresholds but remains short of full acceptance. This reinforces the need to separate optimization-grid dependence from fixed-design solution convergence. The final automated run reports **104 tests passed**.

| Phase 1 requirement | Evidence |
|---|---|
| Numerical implementation runs | Three complete volume-setting bundles with optimizer histories, Q4/P2 fields, diagnostics and cost ledgers |
| Automated verification | [Final test record](../results/tests_final.xml) and [NUMERICAL_VALIDATION.md](NUMERICAL_VALIDATION.md) |
| Reproducible benchmark | Two canonical runs with identical raw/physical arrays and timing-independent numerical history in [summary.json](../results/phase1_v040_verified/summary.json) |
| Verification boundary documented | Nominal numerical acceptance distinguished from angular/stress diagnostics, mass, geometry and physical validation |
| Rerun commands and dependencies | [Root README](../README.md), exact tested versions and immutable manifests/snapshots |
| Research decision stated | Phase 2 remains deferred; the next step is a numerical diagnosis, not an unproven selective algorithm |

The six requested reports are [implementation](IMPLEMENTATION_REPORT.md), [validation](NUMERICAL_VALIDATION.md), [protocol](EXPERIMENT_PROTOCOL.md), [results](BASELINE_RESULTS.md), [limitations](KNOWN_LIMITATIONS.md) and this decision. The Phase 0 package and its references are preserved. The initial Windows telemetry abort and all geometry/resolution failures remain in the evidence.

## Next research gate, not work started in this phase

1. Obtain trustworthy reference labels for the frozen binary candidates. Refine selectively near pixel corners/ligaments or run a separately budgeted finer study; retain the original h2/1/.5 failure statuses. If the geometry changes, assign a new identity and repeat mass, connectivity and complete verification.
2. Diagnose density-to-binary conversion and the 35% disconnected components. Any repair must be a new explicitly costed candidate; do not alter archived failures. Resolve the distinction between a numerical density model and an actual geometry.
3. Before an engineering application, justify the load, material, supports, displacement/strength requirements and model fidelity. The current large-displacement diagnostic is a reason for model review, not an invitation to infer a safe load by extrapolation.
4. Reconsider selective verification only with reliable fine labels, a larger within-constraint catalogue, held-out families and measured cost headroom. Freeze the exact proposed rule and compare it against coarse ordering, seeded random/space-filling, expected verified improvement per cost and the closest finite-catalogue knowledge-gradient baseline. Charge common generation, fitting, failed attempts and mandatory final verification equally.

If conventional geometry checks and deterministic refinement explain the discrepancies or remain cheap enough, keep that workflow and stop the learned-method branch. No Phase 2 implementation, AI integration, UI or infrastructure has begun.
